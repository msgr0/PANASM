import gfapy as gp
import sys
import argparse as ap
from itertools import combinations


"""
this python script removes nodes from a Gfa files based on a provided length threshold
and then reconnects neighbours of the removed node by pairs
"""


def main(args):
    gfa_file_path = args.input
    graph = gp.Gfa.from_file(gfa_file_path)
    remove(graph, args.threshold)
    graph.to_file(args.output)


def remove(gfa, threshold):
    total_len = len(gfa.segments)
    counter = 0
    for seg in gfa.segments:
        counter += 1
        if counter % 10 == 0:
            print(
                "working...",
                int(counter / total_len * 100),
                "% graph analized",
                end="\r",
            )
        # print(seg.name)
        if (seg.LN != None and seg.LN <= threshold) or (len(seg.sequence) <= threshold):
            nodes_to_reconnect = set()

            for e in seg.dovetails:
                if e.from_segment.name != seg.name:
                    nodes_to_reconnect.add(f"{e.from_segment.name}\t{e.from_orient}")
                elif e.to_segment.name != seg.name:
                    nodes_to_reconnect.add(f"{e.to_segment.name}\t{e.to_orient}")
                else:
                    continue  ## self edge
                gfa.rm(e)

            pairs = list(combinations(nodes_to_reconnect, 2))
            for edge in pairs:
                new_edge = gp.Line(f"L\t{edge[0]}\t{edge[1]}\t0M")
                try:
                    gfa.add_line(new_edge)
                except:
                    pass
                    # print(f"edge already added!")
                # print(new_edge)
            seg.disconnect()
            gfa.validate()


if __name__ == "__main__":
    parser = ap.ArgumentParser(
        description="remove gfa segments below the threshold and connect neighbours"
    )
    parser.add_argument("-i", "--input", help="input graph to trim")
    parser.add_argument("-o", "--output", help="output graph")
    parser.add_argument(
        "-t", "--threshold", help="threshold to remove contigs", type=int
    )
    parser.add_argument("-c", "--contains", help="string to remove contigs")

    args = parser.parse_args()

    main(args)
