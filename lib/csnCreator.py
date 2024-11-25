import argparse

parser = argparse.ArgumentParser(description="CSN Creator Script")
parser.add_argument("-m", "--matrix", help="Matrix Path", type=str, required=True)
parser.add_argument(
    "-t", "--threshold", help="CSN Threshold", type=float, required=True
)
parser.add_argument(
    "-o", "--output", help="Output file for the CSN.", type=str, required=True
)

args = parser.parse_args()
_MATRIX = args.matrix
_THRESHOLD = args.threshold
_OUTPUT = args.output


########################################

########################################

import networkx as nx
import pandas as pd


def createCSN():
    mat = pd.read_csv(_MATRIX)
    mat.pop("Unnamed: 0")

    csn = nx.Graph()

    cols = mat.columns
    row_count = 0

    col_arr = []

    for row in mat.itertuples():
        col_count = row_count

        for col in row[row_count + 1 :]:
            if float(col) >= _THRESHOLD:  # If above threshold
                csn.add_edge(cols[row_count], cols[col_count])

                if cols[row_count] not in col_arr:  # For edgeless nodes
                    col_arr += [cols[row_count]]

            col_count += 1

        row_count += 1

    # Add edgeless nodes
    for col in cols:
        if col not in col_arr:
            csn.add_node(col)

    nx.write_graphml(csn, _OUTPUT)


########################################

########################################


def main():
    createCSN()


if __name__ == "__main__":
    main()
