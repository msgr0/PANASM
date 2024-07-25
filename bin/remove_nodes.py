import gfapy as gp
import sys
import argparse as ap
from itertools import combinations


true = True
false = False
null = None


"""
this python script removes nodes from a Gfa files based on a provided length threshold
and then reconnects neighbours of the removed node by pairs
"""


def main(args):
    gfa_file_path = args.input
    graph = gp.Gfa.from_file(gfa_file_path)
    remove(graph, args.threshold)
    graph.to_file(args.output)


def edge_exists(edge, collection):
    for e in collection:
        if edge_compare(edge, e) or edge_compare(edge_mirror(edge), e):
            return true
    return false


def edge_compare(edge1, edge2):
    if str(edge1[0][0]) == str(edge2[0][0]) and str(edge1[1][0]) == str(edge2[1][0]):
        if edge1[0][1] == edge2[0][1] and edge1[1][1] == edge2[1][1]:
            return true
    return false


def edge_mirror(edge):
    """
    from an edge( L(node1, orient1, position1), R(node2, orient2, position2))
    output the mirrored edge (R, L) with reversed orientation
    """
    assert edge[0][2] == "l", "edge[0][2] should be left, 'l'"
    assert edge[1][2] == "r", "edge[1][2] should be right, 'r'"

    return (
        (edge[1][0], invert(edge[1][1]), "l"),
        (edge[0][0], invert(edge[0][1]), "r"),
    )


def edge_tos(edge):
    """
    from an edge( L(node1, orient1, position1), R(node2, orient2, position2))
    output the corresponding gfa formatted Link (L) string.
    """
    # edge[0] L_NODE
    # edge[0][0] name
    # edge[0][1] orient
    # edge[0][2] original pos (l=left, r=right)
    #
    # edge[1] R_NODE
    # edge[1][0] name
    # edge[1][1] orient
    # edge[1][2] original pos (l=left, r=right)
    #
    # if L_NODE(edge[0][2] == 'r') then reverse orient
    # if R_NODE(edge[1][2] == 'l') then reverse orient
    #
    # after that, if edge[0][1] == '-' == edge[1][1] then reverse polarity and invert position
    # i.e.    A -  B -   becomes   B +  A +    for clarity and consistency ...
    #
    l_orient = "?"
    r_orient = "?"

    if edge[0][2] == "r":
        l_orient = invert(edge[0][1])
    elif edge[0][2] == "l":
        l_orient = edge[0][1]

    if edge[1][2] == "l":
        r_orient = invert(edge[1][1])
    elif edge[1][2] == "r":
        r_orient = edge[1][1]

    return gp.Line(f"L\t{edge[0][0]}\t{l_orient}\t{edge[1][0]}\t{r_orient}\t0M")


def invert(sign):
    if sign == "+":
        return "-"
    elif sign == "-":
        return "+"
    if sign == "l":
        return "r"
    elif sign == "r":
        return "l"


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
        #
        edge_collection = set()

        if (seg.LN != None and seg.LN <= threshold) or (len(seg.sequence) <= threshold):
            nodes_to_reconnect = set()

            for e in seg.edges:
                if e.from_segment.name != seg.name:
                    nodes_to_reconnect.add(
                        (f"{e.from_segment.name}", f"{e.from_orient}", "l")
                    )
                elif e.to_segment.name != seg.name:
                    nodes_to_reconnect.add(
                        (f"{e.to_segment.name}", f"{e.to_orient}", "r")
                    )
                gfa.rm(e)  ## else do nothing and remove self edge.

            pairs = list(combinations(nodes_to_reconnect, 2))
            for edge in pairs:
                new_edge = edge_tos(edge)
                # print (new_edge)
                try:
                    gfa.add_line(new_edge)
                except:
                    pass
                    # print(f"edge already added!")
                # print(new_edge)
            gfa.validate()
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
