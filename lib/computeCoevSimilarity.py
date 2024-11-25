import argparse
import networkx as nx
import numpy as np
import pandas as pd
import polars as pl
from itertools import combinations

pl.Config.set_tbl_rows(100)


def compute_score_matrix(coev_path):
    """Computes ECC score matrix prior to jaccard similarity calculation"""

    # Determine cliques in coev graph
    coev_net = nx.read_graphml(coev_path)
    cliques = list(nx.find_cliques(coev_net))
    cliques_int = [list(map(int, clique)) for clique in cliques]

    # Import pair df
    df = pl.read_csv("result_df.csv")

    # Initialise a labelled 2D score matrix
    unique_ids = df.select(pl.col("ID").unique()).to_series().to_list()
    matrix_size = len(unique_ids)
    init_matrix = np.zeros((matrix_size, matrix_size), dtype=int)
    score_matrix = pl.DataFrame(init_matrix, unique_ids)
    score_matrix = score_matrix.with_columns(pl.Series("id", unique_ids))
    score_matrix = score_matrix.select(["id"] + score_matrix.columns[:-1])

    for clique in cliques_int:
        """
        Filter based on current clique coev residues, then group to single
        protein per row
        """
        potential_ecc = df.filter(
            (pl.col("sink_aln").is_in(clique)) | (pl.col("source_aln").is_in(clique))
        )
        grouped = potential_ecc.group_by("ID").agg(
            pl.col("source_aln"), pl.col("sink_aln")
        )

        """
        Determine all-or-nothing Equivalent Coevolving Cliques (ECC). If target
        protein matches exact residue coev pairs then +1 in score matrix, else 
        no increment. Does not allow partial matches
        """
        for target_row in grouped.iter_rows(named=True):
            target_set = set(zip(target_row["source_aln"], target_row["sink_aln"]))
            for comp_row in grouped.iter_rows(named=True):
                comp_set = set(zip(comp_row["source_aln"], comp_row["sink_aln"]))
                if target_row["ID"] != comp_row["ID"]:
                    if target_set == comp_set:
                        score_matrix = score_matrix.with_columns(
                            pl.when(pl.col("id") == comp_row["ID"])
                            .then(pl.col(target_row["ID"]) + 1)
                            .otherwise(pl.col(target_row["ID"]))
                            .alias(target_row["ID"])
                        )
    return score_matrix


def calculate_jaccard(score_matrix):
    """Calculates the jaccard similarity score"""

    unique_ids = score_matrix.select(pl.col("id").unique()).to_series().to_list()

    # Initialise another empty matrix for Jaccard
    matrix_size = len(unique_ids)
    init_matrix = np.zeros((matrix_size, matrix_size), dtype=int)
    jaccard_matrix = pl.DataFrame(init_matrix, unique_ids)
    jaccard_matrix = jaccard_matrix.with_columns(pl.Series("id", unique_ids))
    jaccard_matrix = jaccard_matrix.select(["id"] + jaccard_matrix.columns[:-1])

    # get maximum number of unique combinations
    comparisons = combinations(unique_ids, 2)

    for id_a, id_b in comparisons:
        if id_a == id_b:
            pass

        # Filter to target columns
        targets = score_matrix.select(pl.col(id_a), pl.col(id_b))

        # Calculate intersect
        intersect_list = []
        for intersection in targets.iter_rows():
            intersect_list.append(min(intersection))
        intersect = sum(intersect_list)

        # Calculate set size
        source = score_matrix.select(pl.col(id_a)).sum().item()
        target = score_matrix.select(pl.col(id_b)).sum().item()

        # Calculate union and jaccard
        union = source + target - intersect
        jaccard = intersect / union if union != 0 else 0

        # Update jaccard symetrically
        jaccard_matrix = jaccard_matrix.with_columns(
            pl.when(pl.col("id") == id_a)
            .then(jaccard)
            .otherwise(pl.col(id_b))
            .alias(id_b)
        )
        jaccard_matrix = jaccard_matrix.with_columns(
            pl.when(pl.col("id") == id_b)
            .then(jaccard)
            .otherwise(pl.col(id_a))
            .alias(id_a)
        )

    return jaccard_matrix


def create_csn(jaccard_matrix, threshold):
    G = nx.Graph()

    for col in jaccard_matrix.columns[1:]:
        scores = jaccard_matrix[col].to_list()
        row_ids = jaccard_matrix["id"].to_list()
        G.add_node(col)

        for score, row_id in zip(scores, row_ids):
            if score != 0.0 and score > threshold:
                G.add_edge(col, row_id, similarity=score)

    nx.write_graphml(G, f"csn_{threshold*100}.graphml")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Coev Similarity Compute Script")
    parser.add_argument(
        "-c", "--coev", help="Coevolutionary Graph", type=str, required=True
    )
    parser.add_argument(
        "-t",
        "--threshold",
        help="Threshold for similarity (i.e 0.4 for 40% or more)",
        type=float,
        required=True,
    )
    args = parser.parse_args()

    coev_path = args.coev
    threshold = args.threshold

    score_matrix = compute_score_matrix(coev_path)
    score_matrix.write_csv("score_matrix.csv")

    jaccard = calculate_jaccard(score_matrix)
    jaccard.write_csv("jaccard.csv")

    create_csn(jaccard, threshold)
