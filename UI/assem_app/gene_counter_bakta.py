#!/bin/env python
import os 
import sys
import pandas as pd
import argparse



def gene_counter(gene_files: str, gene_type: str, out_name: str) -> None:

 
    header = ["Sequence_Id",    "Type",    "Start",   "Stop",    "Strand", "Locus_tag",       "Gene",    "Product",    "DbXrefs"]

    gene_counts = {}
    gnt = gene_type.split(',') #to allow multiple genes to be passed on comma-seperated

    for i in os.listdir(gene_files):
        if i.endswith(".tsv"):
            file_path = os.path.join(gene_files, i)
            df = pd.read_csv(file_path, sep= '\t', header=None, names = header, skiprows=1)
            df =df[df['Type'].isin(gnt)]
            name = os.path.splitext(i)[0]
            gene_count = df[["Product"]].value_counts()
            gene_counts[name] = gene_count

    df = pd.DataFrame(gene_counts).fillna(0).astype(int)

    out_path = os.path.join(os.getcwd(), f"{out_name}.tsv")
    df.to_csv(out_path, sep='\t', index=1)
    print(f"Gene count table saved to {out_path}")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog = "Gene counter",
        description= "Takes a list of tsv file of annotated genes and quantify them per genome"
    )


    parser.add_argument("-d", '--gene_files', required=True, help="Path to a directory with tsv files containing gene annotations")
    parser.add_argument('-t', '--gene_type', required=True, help="A list of comma sperated gene types to be counted, e.g cds, tRNA etc.")
    parser.add_argument("-o", '--out_name', default="gene_counts", help="Name of output file. It will be created at the current working directory!")

    args = parser.parse_args()
    gene_counter(args.gene_files, args.gene_type, args.out_name)



