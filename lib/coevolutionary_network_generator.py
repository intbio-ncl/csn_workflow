import numpy as np
import networkx as nx
from copy import deepcopy
from datetime import datetime
from pathlib import Path


def create_network(data_path, file, output_name, node_number):
    """Creates an alignment network of resiude-residue coeveolution positions"""

    temp_mat = np.loadtxt(file)

    indices = [[], []]

    for i in range(0, node_number):  # Picking the top N pairs
        index = np.where(temp_mat == np.max(temp_mat))
        indices[0] += [index[0][0]]
        indices[1] += [index[1][0]]

        temp_mat[index[0][0], index[1][0]] = 0
        temp_mat[index[1][0], index[0][0]] = 0

    zip_obj = zip(indices[0], indices[1])  # Reordering them
    indices[0] = [sorted(point)[0] for point in deepcopy(zip_obj)]
    indices[1] = [sorted(point)[1] for point in zip_obj]

    G = nx.Graph()

    for i in range(len(indices[0])):
        G.add_edge(
            str(indices[0][i] + 1),
            str(indices[1][i] + 1),
        )

    nx.write_graphml(G, f"{str(data_path)}/{output_name}.graphml")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Coevolution Network Creator Script")
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

    parser.add_argument(
        "-o",
        "--output_name",
        help="Output file name of Coevolution graph",
        type=str,
        required=True,
    )
    data_path = Path("./data") / Path(datetime.now().strftime("%y%m%d&m"))

    args = parser.parse_args()
    ccmpred_file = args.file
    coev_graph_name = args.output_name
    node_number = args.num

    create_network(data_path, ccmpred_file, coev_graph_name, node_number)
