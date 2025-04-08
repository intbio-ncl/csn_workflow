import argparse
import multiprocessing
import networkx as nx
import numpy as np
import polars as pl
import time
from multiprocessing import Pool, set_start_method, Manager
import os
from pathlib import Path
from datetime import datetime

pl.Config.set_tbl_rows(100)


def worker(clique, df, matrix_size, unique_ids):
    """Worker abstraction to allow multiprocessing"""
    return clique_scoring(clique, df, matrix_size, unique_ids)


def compute_score_matrix(data_path, coev_path, cpu_n):
    """Computes ECC score matrix prior to jaccard similarity calculation"""
    # Determine cliques in coev graph
    coev_net = nx.read_gml(coev_path)
    cliques = list(nx.find_cliques(coev_net))
    cliques_int = [list(map(int, clique)) for clique in cliques]

    # Import pair df
    df = pl.read_csv(f"{data_path}/result_df.csv")

    # Initialise a labelled 2D score matrix
    unique_ids = df.select(pl.col("ID").unique()).to_series().to_list()
    unique_ids = sorted(unique_ids, key=str.casefold)
    matrix_size = len(unique_ids)

    with Manager() as manager:
        # Prepare inputs for the worker
        inputs = [(clique, df, matrix_size, unique_ids) for clique in cliques_int]

        """ 
        Specifically uses spawn for multiprocessing. Polars is already 
        multithreaded, so forking or alternative multiprocessing techniques
        will not work due to shared state. Neither does locking files due to I/O.
        Spawn, whilst not optimal still leads to significant speed up
        """
        set_start_method("spawn", force=True)
        with Pool(cpu_n) as pool:
            local_matrices = pool.starmap(worker, inputs)
        # Prepare global score matrix
        init_matrix = np.zeros((matrix_size, matrix_size), dtype=int)
        global_matrix = pl.DataFrame(init_matrix, unique_ids)
        global_matrix = global_matrix.with_columns(pl.Series("id", unique_ids))
        global_matrix = global_matrix.select(["id"] + global_matrix.columns[:-1])

        # Merge all score matricies
        for local_matrix in local_matrices:
            global_matrix = merge_matrices(global_matrix, local_matrix)

        global_matrix.write_parquet(f"{data_path}/test.parquet")

    return global_matrix


def clique_scoring(clique, pair_df, matrix_size, unique_ids):
    """Generates scores for all-or-othing Equivalent Coevolving Cliques (ECCs)"""
    now = time.time()

    init_matrix = np.zeros((matrix_size, matrix_size), dtype=int)
    score_matrix = pl.DataFrame(init_matrix, unique_ids)
    score_matrix = score_matrix.with_columns(pl.Series("id", unique_ids))
    score_matrix = score_matrix.select(["id"] + score_matrix.columns[:-1])
    """
    Filter based on current clique coev residues, then group to single
    protein per row. Some Cliques may have been removed from further 
    processing due to their frequency.
    """
    potential_ecc = pair_df.filter(
        (pl.col("sink_aln").is_in(clique)) | (pl.col("source_aln").is_in(clique))
    )
    if len(potential_ecc) == 0:
        pass
    else:
        grouped = potential_ecc.group_by("ID").agg(
            pl.col("source_aln"), pl.col("sink_aln")
        )

        # Combine into a struct
        grouped = grouped.with_columns(
            pl.struct("source_aln", "sink_aln").alias("source_sink")
        )

        # Compare the struct fields
        comparison = grouped.join(grouped, how="cross", suffix="_other").select(
            [
                pl.col("ID").alias("ID_1"),
                pl.col("source_sink").alias("source_sink_1"),
                pl.col("ID_other").alias("ID_2"),
                pl.col("source_sink_other").alias("source_sink_2"),
                (pl.col("source_sink") == pl.col("source_sink_other"))
                .cast(pl.Int32)
                .alias("binary_score"),
            ]
        )

        unique_combinations = comparison.filter(
            (pl.col("ID_1") < pl.col("ID_2")) & (pl.col("ID_1") != pl.col("ID_2"))
        )

        scores = unique_combinations.filter(pl.col("binary_score") == 1)
        """
        Determine all-or-nothing Equivalent Coevolving Cliques (ECC). If target
        protein matches exact residue coev pairs then +1 in score matrix, else 
        no increment. Does not allow partial matches.
        """
        # Extract indices for ID_1 and ID_2
        for row in scores.filter(pl.col("binary_score") == 1).iter_rows(named=True):
            id_1 = row["ID_1"]
            id_2 = row["ID_2"]
            score_matrix = score_matrix.with_columns(
                pl.when(pl.col("id") == id_1)
                .then(pl.col(id_2) + 1)
                .otherwise(pl.col(id_2))
                .alias(id_2)
            )

            score_matrix = score_matrix.with_columns(
                pl.when(pl.col("id") == id_2)
                .then(pl.col(id_1) + 1)
                .otherwise(pl.col(id_1))
                .alias(id_1)
            )

    time_taken = time.time() - now

    print(f" Time: {time_taken:.3f} seconds \t Clique: {clique}")

    return score_matrix


def merge_matrices(global_matrix, local_matrix):
    """Merge a thread-local score matrix into the global matrix"""

    for col in local_matrix.columns[1:]:
        global_matrix = global_matrix.with_columns(
            (pl.col(col) + local_matrix[col]).alias(col)
        )

    return global_matrix


def calculate_jaccard(score_matrix):
    """Calculates Jaccard similarity"""

    # Extract IDs and convert the numeric data to a NumPy array
    ids = score_matrix["id"].to_list()
    matrix = score_matrix.drop("id").to_numpy()

    # Compute the intersection for all pairs: element-wise min across rows
    intersection = np.minimum(matrix[:, :, None], matrix[:, None, :]).sum(axis=0)

    # Compute the union for all pairs: row sums minus the intersection
    row_sums = matrix.sum(axis=0)
    union = row_sums[:, None] + row_sums[None, :] - intersection

    # Compute Jaccard index, avoiding division by zero
    jaccard_matrix = np.divide(
        intersection,
        union,
        out=np.zeros_like(intersection, dtype=float),
        where=union != 0,
    )

    jaccard_df = (
        pl.DataFrame(jaccard_matrix, schema=ids)
        .with_columns(pl.Series("id", ids))
        .melt(id_vars="id", variable_name="id2", value_name="jaccard")
    )

    print("Jaccard similarity computed")
    return jaccard_df


def create_csn(data_path, jaccard_matrix, threshold):
    """Creates Coevolution Similarity Network (CSN) with connected nodes above threshold"""

    print("Creating CSN")
    G = nx.Graph()

    source = jaccard_matrix["id"].to_list()
    target = jaccard_matrix["id2"].to_list()
    score = jaccard_matrix["jaccard"].to_list()

    for target_threshold in threshold:
        for v1, v2, e in zip(source, target, score):
            G.add_node(v1)
            if e != 0.0 and e > target_threshold:
                G.add_edge(v1, v2, similarity=e)

        nx.write_gml(G, f"{data_path}/csn_vec_{target_threshold * 100}.gml")
        print("CSN has been created")


def compute_coevolutionary_similarity(data_path, coev_path, threshold, cpu_n):
    if cpu_n > multiprocessing.cpu_count():
        raise RuntimeError("cpu count greater than physical cores available")

    if not os.path.exists(f"{data_path}/test.parquet"):
        score_matrix = compute_score_matrix(data_path, coev_path, cpu_n)
        score_matrix.write_csv(f"{data_path}/score.csv")
    else:
        score_matrix = pl.read_parquet(f"{data_path}/test.parquet")

    jaccard = calculate_jaccard(score_matrix)
    jaccard.write_csv(f"{data_path}/jaccard.csv")

    create_csn(data_path, jaccard, threshold)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Coev Similarity Compute Script")
    parser.add_argument(
        "-cg", "--coev", help="Coevolutionary Graph", type=str, required=True
    )
    parser.add_argument(
        "-c",
        "--cpu",
        help="Number of cores to use for multiprocessing",
        type=int,
        required=True,
    )
    args = parser.parse_args()
    data_path = Path.cwd() / Path("test")

    coev_path = args.coev
    cpu_n = args.cpu
    if cpu_n > multiprocessing.cpu_count():
        raise RuntimeError("cpu count greater than physical cores available")

    if not os.path.exists(f"{data_path}/test.parquet"):
        score_matrix = compute_score_matrix(data_path, coev_path, cpu_n)
        score_matrix.write_csv(f"{data_path}/score.csv")
    else:
        score_matrix = pl.read_parquet(f"{data_path}/test.parquet")

    jaccard = calculate_jaccard(score_matrix)
    jaccard.write_csv(f"{data_path}/jaccard.csv")

    threshold = [0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.975, 0.99, 0.995]
    create_csn(data_path, jaccard, threshold)
