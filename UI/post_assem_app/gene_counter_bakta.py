#!/bin/env python
import os
import sys
import pandas as pd
import argparse


def _annotation_tables(gene_files: str):
    path = os.path.abspath(gene_files)
    tables = []
    if os.path.isfile(path):
        parent = os.path.dirname(path)
        if path.lower().endswith((".tsv", ".txt")):
            sibling_dir = parent if os.path.isdir(parent) else None
        else:
            sibling_dir = None
        search_dir = sibling_dir
        if search_dir:
            for name in os.listdir(search_dir):
                if name.lower().endswith((".tsv", ".txt")):
                    tables.append(os.path.join(search_dir, name))
        if not tables:
            tables.append(path)
    elif os.path.isdir(path):
        for root, _dirs, names in os.walk(path):
            for name in names:
                if name.lower().endswith((".tsv", ".txt")):
                    tables.append(os.path.join(root, name))
    return sorted(set(tables))


def gene_counter(gene_files: str, gene_type: str, out_name: str) -> None:
    header = ["Sequence_Id", "Type", "Start", "Stop", "Strand", "Locus_tag", "Gene", "Product", "DbXrefs"]
    gene_counts = {}
    gnt = [g.strip() for g in str(gene_type).split(",") if g.strip()]
    gnt_l = {g.lower() for g in gnt}

    tables = _annotation_tables(gene_files)
    if not tables:
        raise SystemExit(f"No .tsv/.txt annotation tables found at {gene_files}")

    for file_path in tables:
        df = pd.read_csv(file_path, sep="\t", header=None, names=header, skiprows=1)
        if "Type" in df.columns:
            df = df[df["Type"].astype(str).str.lower().isin(gnt_l)]
        name = os.path.splitext(os.path.basename(file_path))[0]
        gene_count = df[["Product"]].value_counts()
        gene_counts[name] = gene_count

    df = pd.DataFrame(gene_counts).fillna(0).astype(int)
    out_path = f"{out_name}.tsv" if not str(out_name).endswith(".tsv") else out_name
    if not os.path.isabs(out_path) and os.sep not in str(out_name) and os.sep not in out_path:
        out_path = os.path.join(os.getcwd(), out_path)
    parent = os.path.dirname(out_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=1)
    print(f"Gene count table saved to {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Gene counter",
        description="Takes a list of tsv file of annotated genes and quantify them per genome",
    )
    parser.add_argument("-d", "--gene_files", required=True, help="Path to a directory or a single TSV/TXT annotation file")
    parser.add_argument("-t", "--gene_type", required=True, help="Comma-separated gene types to count, e.g. cds,tRNA")
    parser.add_argument("-o", "--out_name", default="gene_counts", help="Output file name or path")
    args = parser.parse_args()
    gene_counter(args.gene_files, args.gene_type, args.out_name)
