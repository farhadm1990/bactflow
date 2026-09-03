#!/bin/env python
import os
import sys
import pandas as pd
import argparse


BAKTA_COLS = [
    "Sequence_Id",
    "Type",
    "Start",
    "Stop",
    "Strand",
    "Locus_tag",
    "Gene",
    "Product",
    "DbXrefs",
]


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


def _read_bakta_table(file_path: str) -> pd.DataFrame:
    """Parse Bakta TSV: skip # comments, use commented header when present."""
    header = None
    rows = []
    with open(file_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                body = line[1:].strip()
                # Bakta header line, e.g. "#Sequence Id\tType\t...\tGene\tProduct\tDbXrefs"
                if "\t" in body and "Product" in body and "Type" in body:
                    header = [h.strip().replace(" ", "_") for h in body.split("\t")]
                continue
            if not line.strip():
                continue
            rows.append(line.rstrip("\n").split("\t"))

    if not rows:
        return pd.DataFrame(columns=BAKTA_COLS)

    width = max(len(r) for r in rows)
    if header and len(header) >= 8:
        cols = header[:width]
        if len(cols) < width:
            cols = cols + [f"col_{i}" for i in range(len(cols), width)]
    else:
        cols = BAKTA_COLS[:width] if width <= len(BAKTA_COLS) else BAKTA_COLS + [
            f"col_{i}" for i in range(len(BAKTA_COLS), width)
        ]

    padded = [r + [""] * (len(cols) - len(r)) for r in rows]
    df = pd.DataFrame(padded, columns=cols)

    # Normalize expected names if Bakta used spaces
    rename = {}
    for col in df.columns:
        key = col.replace(" ", "_")
        if key != col:
            rename[col] = key
    if rename:
        df = df.rename(columns=rename)

    for needed in BAKTA_COLS:
        if needed not in df.columns:
            # try loose match
            lower = {c.lower(): c for c in df.columns}
            if needed.lower() in lower:
                df = df.rename(columns={lower[needed.lower()]: needed})
            elif needed == "Locus_tag" and "Locus_Tag" in df.columns:
                df = df.rename(columns={"Locus_Tag": "Locus_tag"})
            else:
                df[needed] = ""

    return df[BAKTA_COLS]


def gene_counter(gene_files: str, gene_type: str, out_name: str) -> None:
    gnt = [g.strip() for g in str(gene_type).split(",") if g.strip()]
    gnt_l = {g.lower() for g in gnt}

    tables = _annotation_tables(gene_files)
    if not tables:
        raise SystemExit(f"No .tsv/.txt annotation tables found at {gene_files}")

    gene_counts = {}
    map_frames = []

    for file_path in tables:
        df = _read_bakta_table(file_path)
        if gnt_l:
            df = df[df["Type"].astype(str).str.lower().isin(gnt_l)]
        name = os.path.splitext(os.path.basename(file_path))[0]
        products = df["Product"].fillna("").astype(str).str.strip()
        products = products[products != ""]
        gene_counts[name] = products.value_counts()

        genes = df["Gene"].fillna("").astype(str).str.strip()
        map_frames.append(
            pd.DataFrame(
                {
                    "Genome": name,
                    "Gene": genes,
                    "Product": df["Product"].fillna("").astype(str).str.strip(),
                }
            )
        )

    count_df = pd.DataFrame(gene_counts).fillna(0).astype(int)
    count_df.index.name = "Product"

    out_path = f"{out_name}.tsv" if not str(out_name).endswith(".tsv") else out_name
    if not os.path.isabs(out_path) and os.sep not in str(out_name) and os.sep not in out_path:
        out_path = os.path.join(os.getcwd(), out_path)
    parent = os.path.dirname(out_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    count_df.to_csv(out_path, sep="\t", index=True)
    print(f"Gene count table saved to {out_path}")

    # Gene↔Product map so the plotter can match enzyme names against Bakta Gene and Product
    map_path = out_path[:-4] + ".gene_map.tsv" if out_path.endswith(".tsv") else out_path + ".gene_map.tsv"
    map_df = pd.concat(map_frames, ignore_index=True)
    map_df = map_df[(map_df["Product"] != "") | (map_df["Gene"] != "")]
    map_df = map_df.drop_duplicates()
    map_df.to_csv(map_path, sep="\t", index=False)
    print(f"Gene–product map saved to {map_path}")


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
