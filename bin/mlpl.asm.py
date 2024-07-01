import argparse
import gfapy as gp
import numpy
import pandas as pd
import sys


def main(args):
    pred = args.pred
    graph = args.graph
    output = args.output
    pbf = args.pbf
    prefix = args.prefix

    mlpred_head = [
        "Prob_Chromosome",
        "Prob_Plasmid",
        "Prediction",
        "Contig_name",
        "Contig_length",
    ]

    mlpred = pd.read_csv(pred, sep="\t", header=0)
    mlpred["Contig_name"] = mlpred["Contig_name"].astype(str)
    # return
    df = pd.DataFrame(columns=mlpred_head)

    df_pbf = pd.DataFrame(columns=["contig", "pscore"])

    gfa = gp.Gfa.from_file(graph)

    for frag in gfa.segments:
        # contig_list = frag.cl.split(",")
        contig = frag.name  # DOESNT HAVE uni or ske if GPLAS
        # ncontigs = len(contig_list)
        prob_chromo = 0
        prob_plas = 0
        length = frag.LN if frag.LN != None else len(frag.sequence)
        # for contig in contig_list:
        # print(contig)
        # eliminated ske-uni prefix
        try:
            mlpred.loc[
                mlpred["Contig_name"] == f"{prefix}{contig}", "Prediction"
            ].values[0]

        except Exception as e:
            continue

        label = mlpred.loc[
            mlpred["Contig_name"] == str(f"{prefix}{contig}"), "Prediction"
        ].values[0]
        prob_chromo = mlpred.loc[
            mlpred["Contig_name"] == str(f"{prefix}{contig}"), "Prob_Chromosome"
        ].values[0]
        prob_plas = mlpred.loc[
            mlpred["Contig_name"] == str(f"{prefix}{contig}"), "Prob_Plasmid"
        ].values[0]
        # labeling = mlpred.loc[mlpred[3] == contig, 2].values[0]

        name = f"S{contig}_LN:i:{frag.LN}_dp:f:{frag.dp}"

        df.loc[len(df)] = {
            "Prob_Chromosome": prob_chromo,
            "Prob_Plasmid": prob_plas,
            "Prediction": label,
            "Contig_name": name,
            "Contig_length": length,
        }

        if pbf != None:
            df_pbf.loc[len(df_pbf)] = {"contig": contig, "pscore": prob_plas}

    if pbf != None:
        df_pbf.to_csv(pbf, sep="\t", index=False, header=False)
    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="convert mlplasmid pred")

    parser.add_argument("--pred", "-p", help="mlplasmid asm pred")
    parser.add_argument("--graph", "-g", help="pangenome_graph convert mlplas into")
    parser.add_argument("--prefix", help="prefix to add when matching, used for gplas")
    parser.add_argument("--output", "-o", help="gplas output mode")
    parser.add_argument("--pbf", help="pbf output mode")
    args = parser.parse_args()
    main(args)
