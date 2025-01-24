import argparse
from ctypes import alignment
from Bio import AlignIO
import networkx as nx
from itertools import combinations
import polars as pl
from pathlib import Path
from datetime import datetime


def createAlnVec(seq):
    """Vec where index = sequence position and value = alignment position"""

    aln_vec = []
    gapped_seq = str(seq.seq)

    for i in range(len(gapped_seq)):
        if gapped_seq[i] == "-":
            continue

        aln_vec.append(i + 1)

    return aln_vec


def readCoevNetworkVec(seq, G, aln_vec, pair_dict):
    """outputs network as duplex tuple vec, (res_aln, res_prot)(res_aln, res_prot)"""

    source_aln_vec = []
    source_aa_vec = []
    sink_aln_vec = []
    sink_aa_vec = []

    for edge in G.edges():
        source = int(edge[0].split("-")[-1])
        sink = int(edge[1].split("-")[-1])

        # swap so greater residue is always last
        if source > sink:
            temp = sink
            sink = source
            source = temp

        try:
            # Ensures that ValueError is captured prior to any list extension
            source_index = aln_vec.index(source) + 1
            sink_index = aln_vec.index(sink) + 1

            source_aln_vec.append(source)
            source_aa_vec.append(source_index)
            sink_aln_vec.append(sink)
            sink_aa_vec.append(sink_index)
        except ValueError:
            continue

    # Ensures consistency and captures shape error
    assert (
        len(source_aln_vec)
        == len(source_aa_vec)
        == len(sink_aln_vec)
        == len(sink_aa_vec)
    ), "Mismatch in lengths of vectors!"

    pair_dict[seq.id] = {
        "source_aln": source_aln_vec,
        "source_aa": source_aa_vec,
        "sink_aln": sink_aln_vec,
        "sink_aa": sink_aa_vec,
    }

    return pair_dict


def construct_df(data_path, pair_dict):
    """Creates a polars dataframe with ID seq.id and then 4 columns containing
    residue information"""
    frames = []

    for key, nested_data in pair_dict.items():
        df = pl.DataFrame(nested_data)

        df = df.with_columns(pl.lit(key).alias("ID"))

        frames.append(df)

    result = pl.concat(frames, how="vertical")

    result = result.select(["ID"] + [col for col in result.columns if col != "ID"])
    result.write_csv(f"{data_path}/id_df.csv")

    return result


def compute_frequencies(data_path, df):
    """Computes how often coev resiudes occur"""

    frequency_df = df.group_by(["source_aln", "sink_aln"]).agg(
        [
            pl.col("ID").count().alias("frequency"),
        ]
    )

    frequency_df = frequency_df.sort("frequency", descending=True)
    frequency_df.write_csv(f"{data_path}/freq.csv")

    return frequency_df


def remove_noise(data_path, df, cutoff):
    """Removes coevs if they occur in n% of sequences"""

    filtered_df = df.filter(pl.col("frequency") <= cutoff)
    filtered_df.write_csv(f"{data_path}/filtered_df.csv")

    return filtered_df


def createALNGraphDf(data_path, freq_df, full_df):
    """Creates alignment network from dataframe"""

    G = nx.Graph()

    # Perform an inner join on source_aa and sink_aa to find matches
    matched_df = freq_df.join(full_df, on=["source_aln", "sink_aln"], how="inner")

    # Create tuples of source_aln and sink_aln
    result_df = matched_df.select(
        [
            pl.col("ID"),
            pl.col("source_aln"),
            pl.col("sink_aln"),
            pl.col("source_aa"),
            pl.col("sink_aa"),
        ]
    )
    result_df.write_csv(f"{data_path}/result_df.csv")
    G = extract_links(result_df, G)
    # filtered_rows = result_df.filter(pl.col("sink_aln") == 420)

    return G


def extract_links(df, G):
    """
    Extract links for every unique (source_aln, sink_aln) pair.

    """
    # Find unique (source_aln, sink_aln) pairs
    unique_pairs = (
        df.select(["source_aln", "sink_aln"])
        .unique()
        .sort(by="source_aln", descending=False)
    )
    # Loop through each unique (source_aln, sink_aln) pair
    for pair in unique_pairs.iter_rows():
        source_aln, sink_aln = pair

        # Filter rows matching the current pair
        matching_rows = df.filter(
            (pl.col("source_aln") == source_aln) & (pl.col("sink_aln") == sink_aln)
        )

        # Extract `source_aa` and `sink_aa` with `ID` prepended
        source_aa_list = [
            f"{row['ID']}-{row['source_aa']}" for row in matching_rows.to_dicts()
        ]
        sink_aa_list = [
            f"{row['ID']}-{row['sink_aa']}" for row in matching_rows.to_dicts()
        ]

        G.add_edges_from(combinations(source_aa_list, 2))
        G.add_edges_from(combinations(sink_aa_list, 2))

    return G


def create_alignment_network(
    data_path, coev_cutoff, coev_graph_path, aln_graph_name, alignment_file
):
    alignment = AlignIO.read(alignment_file, "clustal")
    seq_number = len(alignment)
    cutoff = seq_number * coev_cutoff
    pair_dict = {}

    G = nx.read_graphml(coev_graph_path)
    for x in range(seq_number):
        current_seq = alignment[x]
        aln_vec = createAlnVec(current_seq)
        pair_vec = readCoevNetworkVec(current_seq, G, aln_vec, pair_dict)

    df = construct_df(data_path, pair_dict)

    frequencies = compute_frequencies(data_path, df)
    cleaned = remove_noise(data_path, frequencies, cutoff)

    G = createALNGraphDf(data_path, cleaned, df)

    print(f"Writing ALN Graph to {aln_graph_name}")
    nx.write_graphml(G, f"{data_path}/{aln_graph_name}.graphml")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Alignment Network Creator Script")
    parser.add_argument(
        "-f",
        "--filter",
        help="Filter % to ignore coevolving columns that occur too frequently (i.e 0.6 = 60%)",
        type=float,
    )
    parser.add_argument("-a", "--aln", help="Alignment file path", type=str)
    parser.add_argument("-g", "--coev", help="Coev graph file path", type=str)
    parser.add_argument("-o", "--output", help="Alignment graph file name", type=str)
    args = parser.parse_args()
    data_path = Path("./data") / Path(datetime.now().strftime("%y%m%d"))

    coev_cutoff = args.filter
    coev_graph_path = args.coev
    aln_graph_name = args.output
    alignment_file = args.aln

    alignment = AlignIO.read(alignment_file, "clustal")
    seq_number = len(alignment)
    cutoff = seq_number * coev_cutoff
    pair_dict = {}

    G = nx.read_graphml(coev_graph_path)
    for x in range(seq_number):
        current_seq = alignment[x]
        aln_vec = createAlnVec(current_seq)
        pair_vec = readCoevNetworkVec(current_seq, G, aln_vec, pair_dict)

    df = construct_df(data_path, pair_dict)

    frequencies = compute_frequencies(data_path, df)
    cleaned = remove_noise(data_path, frequencies, cutoff)

    G = createALNGraphDf(data_path, cleaned, df)

    print(f"Writing ALN Graph to {aln_graph_name}")
    nx.write_graphml(G, f"{data_path}/{aln_graph_name}.graphml")
