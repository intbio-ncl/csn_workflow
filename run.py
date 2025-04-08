from lib.coevolutionary_network_generator import create_network
from datetime import datetime
from lib.alignment_network_generator import create_alignment_network
from lib.coev_similarity_network_generator import compute_coevolutionary_similarity
import argparse
from pathlib import Path
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate co-evolutionary similarity networks"
    )
    parser.add_argument(
        "-n",
        "--num",
        help="Number of nodes for the Coevolution Networks",
        type=int,
        required=True,
    )
    parser.add_argument(
        "-f",
        "--file",
        help="Path to file/folder of the Coevolution Matrices",
        type=str,
        required=True,
    )
    parser.add_argument("-a", "--aln", help="Alignment file path", type=str)
    parser.add_argument(
        "-fi",
        "--filter",
        help="Filter % to ignore coevolving columns that occur too frequently (i.e 0.6 = 60%)",
        type=float,
        default=0.6,
    )
    parser.add_argument(
        "-t",
        "--threshold",
        help="Threshold for similarity (i.e 0.4 for 40% or more)",
        type=float,
        default=0.4,
    )
    parser.add_argument(
        "-c",
        "--cpu",
        help="Number of cores to use for multiprocessing",
        type=int,
        default=4,
    )
    cwd = os.getcwd()
    data_path = Path(datetime.now().strftime("%y%m%d%H"))

    if not data_path.exists():
        # Create the directory (including intermediate directories if needed)
        data_path.mkdir(parents=True)

    args = parser.parse_args()
    ccmpred_file = args.file
    node_number = args.num
    coev_cutoff = args.filter
    alignment_file = args.aln
    threshold = args.threshold
    cpu_n = args.cpu

    create_network(data_path, ccmpred_file, "coevolutionary_network", node_number)
    create_alignment_network(
        data_path,
        coev_cutoff,
        f"{data_path}/coevolutionary_network.graphml",
        "alignment_network",
        alignment_file,
    )
    compute_coevolutionary_similarity(
        data_path, f"{data_path}/coevolutionary_network.graphml", threshold, cpu_n
    )
