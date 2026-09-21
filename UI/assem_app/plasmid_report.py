#!/usr/bin/env python3
"""Load geNomad plasmid summaries for the assembly UI."""

from __future__ import annotations

import csv
import os

DISPLAY_FIELDS = [
    ("Genome", "Genome"),
    ("seq_name", "Contig"),
    ("length", "Length"),
    ("topology", "Topology"),
    ("plasmid_score", "Score"),
    ("n_hallmarks", "Hallmarks"),
    ("conjugation_genes", "Conjugation genes"),
    ("amr_genes", "AMR genes"),
]


def find_plasmid_summary(out_dir):
    if not out_dir:
        return ""
    candidates = [
        os.path.join(out_dir, "plasmid_out", "plasmid_summary.tsv"),
        os.path.join(out_dir, "plasmid_summary.tsv"),
    ]
    for path in candidates:
        if os.path.isfile(path) and os.path.getsize(path) > 0:
            return path
    root = os.path.join(out_dir, "plasmid_out")
    if os.path.isdir(root):
        for dirpath, _dirs, files in os.walk(root):
            for name in files:
                if name.endswith("plasmid_summary.tsv"):
                    path = os.path.join(dirpath, name)
                    if os.path.isfile(path) and os.path.getsize(path) > 0:
                        return path
    return ""


def _pretty(value):
    text = str(value or "").strip()
    if text.lower() in ("", "nan", "none", "na", "."):
        return ""
    return text


def _score(value):
    text = _pretty(value)
    if not text:
        return ""
    try:
        return round(float(text), 3)
    except (TypeError, ValueError):
        return text


def plasmid_display_rows(out_dir):
    path = find_plasmid_summary(out_dir)
    if not path:
        return []
    fallback_genome = os.path.basename(path)
    for suffix in (".plasmid_summary.tsv", "_plasmid_summary.tsv"):
        if fallback_genome.endswith(suffix):
            fallback_genome = fallback_genome[: -len(suffix)]
            break
    rows = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for rec in reader:
            contig = _pretty(rec.get("seq_name") or rec.get("Contig"))
            if not contig:
                continue
            row = {}
            for src, dest in DISPLAY_FIELDS:
                if dest == "Score":
                    row[dest] = _score(rec.get(src) or rec.get("Score"))
                elif dest == "Length":
                    raw = _pretty(rec.get(src) or rec.get("Length"))
                    try:
                        row[dest] = int(float(raw))
                    except (TypeError, ValueError):
                        row[dest] = raw
                elif dest == "Genome":
                    row[dest] = _pretty(rec.get(src) or rec.get(dest)) or fallback_genome
                else:
                    row[dest] = _pretty(rec.get(src) or rec.get(dest))
            rows.append(row)
    return rows
