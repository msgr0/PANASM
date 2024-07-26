"""
Filtering the contigs of length below a spcified threshold in a GFA graph.

Usage:
python filter_GFA.py <input GFA graph> <contig length threshold> <output GFA graph>
"""


import itertools
import networkx as nx
import fileinput
import sys

def get_extremities(v1,s1,v2,s2):
    if s1 == "+":
        w1 = f"{v1}_h" # HEAD of v1 if s1 is plus
    else:
        w1 = f"{v1}_t"
    if s2 == "+":
        w2 = f"{v2}_t"
    else:
        w2 = f"{v2}_h"
    return((w1,w2))

def extract_edge(edge):
    v1,v2 = edge[0],edge[1]
    ext1,ext2 = edge[2][0],edge[2][1]
    if ext1 == f"{v1}_h" or ext2 == f"{v1}_h": s1 = "+"
    elif ext1 == f"{v1}_t" or ext2 == f"{v1}_t": s1 = "-"
    if ext1 == f"{v2}_h" or ext2 == f"{v2}_h": s2 = "-"
    elif ext1 == f"{v2}_t" or ext2 == f"{v2}_t": s2 = "+"
    return((v1,s1,v2,s2))

def filter_graph(in_GFA_file, out_GFA_file, minimum_contig_length):
    """Read a single graph from gfa of gfa.gz, compute attributes, add its nodes and edges to nx graph.
    Contigs shorter than minimum_contig_length are contracted."""

    graph = nx.Graph()

    # first pass: read all nodes
    attributes_nodes = {}
    all_nodes = []

    with fileinput.input(
            in_GFA_file, openhook=fileinput.hook_compressed, mode='r'
    ) as file:
        for line in file:
            if isinstance(line, bytes):
                line = line.decode("utf-8") # convert byte sequences to strings
            parts = line.strip().split("\t")
            if parts[0] == "S":  # node line
                node_id = parts[1]
                seq_len = len(parts[2])

                all_nodes.append(node_id)
                attributes_nodes[node_id] = line.strip()
                assert node_id not in graph
                graph.add_node(node_id)
                graph.nodes[node_id]["length"] = seq_len

    # second pass: read all edges
    with fileinput.input(
            in_GFA_file, openhook=fileinput.hook_compressed, mode='r'
    ) as file:
        for line in file:
            if isinstance(line, bytes):
                line = line.decode("utf-8")
            parts = line.strip().split("\t")
            if parts[0] == "L":  # edge line
                graph.add_edge(
                    parts[1],parts[3],
                    extremities=get_extremities(parts[1],parts[2],parts[3],parts[4])
                )
                # extremities = contigs suffixed by _h for head and _t for tail

    # Filtering short nodes and adding edges betwee extremities of their neighbours
    for node_id in all_nodes:
        if graph.nodes[node_id]["length"] < minimum_contig_length:
            edges = graph.edges(nbunch=node_id, data="extremities")
            neighbours = []
            for edge in edges:
                extremities = edge[2]
                if extremities[0] in [f"{node_id}_h", f"{node_id}_t"]:
                    neighbours.append(extremities[1])
                else:
                    neighbours.append(extremities[0])
            all_new_edges = []
            for (ngb_ext1,ngb_ext2) in list(itertools.combinations(neighbours, 2)):
                graph.add_edge(
                    ngb_ext1[:-2], ngb_ext2[:-2], extremities=(ngb_ext1,ngb_ext2)
                )
            graph.remove_node(node_id)

    # Writing a new GFA file
    with open(out_GFA_file, "w") as out_file:
        for node_id in graph.nodes:
            out_file.write(f"{attributes_nodes[node_id]}\n")
        for edge in graph.edges(data="extremities"):
            print(edge)
            (v1,s1,v2,s2) = extract_edge(edge)
            print (v1, s1, v2, s2)
            edge_str = "\t".join([v1,s1,v2,s2,"0M"])
            out_file.write(f"L\t{edge_str}\n")


filter_graph(sys.argv[1], sys.argv[3], int(sys.argv[2]))
