#!/usr/bin/env python3
"""Load geNomad plasmid + virus summaries for the BactFlow UI (assembly + post-assembly)."""

from __future__ import annotations

import csv
import os
from collections import Counter

PLASMID_DISPLAY_FIELDS = [
    ("Genome", "Genome"),
    ("seq_name", "Contig"),
    ("length", "Length"),
    ("topology", "Topology"),
    ("plasmid_score", "Score"),
    ("n_hallmarks", "Hallmarks"),
    ("conjugation_genes", "Conjugation genes"),
    ("amr_genes", "AMR genes"),
]

VIRUS_DISPLAY_FIELDS = [
    ("Genome", "Genome"),
    ("seq_name", "Contig"),
    ("length", "Length"),
    ("topology", "Topology"),
    ("coordinates", "Coordinates"),
    ("virus_score", "Score"),
    ("n_hallmarks", "Hallmarks"),
    ("marker_enrichment", "Marker enrichment"),
    ("taxonomy", "Taxonomy"),
    ("n_genes", "Genes"),
]


def _find_summary(out_dir, filenames, walk_suffix):
    if not out_dir:
        return ""
    candidates = [os.path.join(out_dir, *parts) for parts in filenames]
    for path in candidates:
        if os.path.isfile(path) and os.path.getsize(path) > 0:
            return path
    root = os.path.join(out_dir, "plasmid_out")
    if os.path.isdir(root):
        for dirpath, _dirs, files in os.walk(root):
            for name in files:
                if name.endswith(walk_suffix):
                    path = os.path.join(dirpath, name)
                    if os.path.isfile(path) and os.path.getsize(path) > 0:
                        return path
    return ""


def find_plasmid_summary(out_dir):
    return _find_summary(
        out_dir,
        [("plasmid_out", "plasmid_summary.tsv"), ("plasmid_summary.tsv",)],
        "plasmid_summary.tsv",
    )


def find_virus_summary(out_dir):
    return _find_summary(
        out_dir,
        [("plasmid_out", "virus_summary.tsv"), ("virus_summary.tsv",)],
        "virus_summary.tsv",
    )


def genomad_scanned_genomes(out_dir):
    """Genome names scanned by geNomad (from genomes_scanned.txt or per-genome summaries)."""
    if not out_dir:
        return []
    names = []
    seen = set()
    for path in (
        os.path.join(out_dir, "plasmid_out", "genomes_scanned.txt"),
        os.path.join(out_dir, "genomes_scanned.txt"),
    ):
        if os.path.isfile(path) and os.path.getsize(path) > 0:
            with open(path, encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    genome = line.strip()
                    if genome and genome not in seen:
                        seen.add(genome)
                        names.append(genome)
            if names:
                return names
    sum_dir = os.path.join(out_dir, "plasmid_out", "summaries")
    if not os.path.isdir(sum_dir):
        sum_dir = os.path.join(out_dir, "summaries")
    if os.path.isdir(sum_dir):
        for name in sorted(os.listdir(sum_dir)):
            genome = ""
            if name.endswith(".plasmid_summary.tsv"):
                genome = name[: -len(".plasmid_summary.tsv")]
            elif name.endswith(".virus_summary.tsv"):
                genome = name[: -len(".virus_summary.tsv")]
            if genome and genome not in seen:
                seen.add(genome)
                names.append(genome)
    if names:
        return names
    for finder in (find_plasmid_summary, find_virus_summary):
        path = finder(out_dir)
        if not path:
            continue
        try:
            with open(path, encoding="utf-8", errors="replace") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for rec in reader:
                    genome = _pretty(rec.get("Genome"))
                    if genome and genome not in seen:
                        seen.add(genome)
                        names.append(genome)
        except OSError:
            continue
    return names


def genomad_run_complete(out_dir):
    """True when geNomad published at least one summary file (even if empty of hits)."""
    return bool(find_plasmid_summary(out_dir) or find_virus_summary(out_dir))


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


def _as_float(value):
    text = _pretty(value)
    if not text:
        return None
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def _as_int(value):
    num = _as_float(value)
    if num is None:
        return None
    try:
        return int(num)
    except (TypeError, ValueError):
        return None


def _display_rows(path, fields, genome_suffixes):
    if not path:
        return []
    fallback_genome = os.path.basename(path)
    for suffix in genome_suffixes:
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
            for src, dest in fields:
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
                elif dest == "Marker enrichment":
                    row[dest] = _score(rec.get(src) or rec.get(dest))
                else:
                    row[dest] = _pretty(rec.get(src) or rec.get(dest))
            rows.append(row)
    return rows


def plasmid_display_rows(out_dir):
    return _display_rows(
        find_plasmid_summary(out_dir),
        PLASMID_DISPLAY_FIELDS,
        (".plasmid_summary.tsv", "_plasmid_summary.tsv"),
    )


def virus_display_rows(out_dir):
    return _display_rows(
        find_virus_summary(out_dir),
        VIRUS_DISPLAY_FIELDS,
        (".virus_summary.tsv", "_virus_summary.tsv"),
    )


def _chart_payload(rows, kind):
    empty = {
        "summary": {
            "n_items": 0,
            "n_genomes": 0,
            "mean_score": None,
            "total_bp": 0,
            "with_amr": 0,
            "with_conjugation": 0,
            "with_taxonomy": 0,
        },
        "by_genome": {"genomes": [], "counts": [], "total_bp": []},
        "scores": {"values": [], "labels": []},
        "lengths": {"contigs": [], "lengths": [], "genomes": []},
        "topology": {"labels": [], "counts": []},
        "taxonomy": {"labels": [], "counts": []},
        "kind": kind,
    }
    if not rows:
        # Keep plasmid-specific keys for older JS
        if kind == "plasmid":
            empty["summary"]["n_plasmids"] = 0
        else:
            empty["summary"]["n_viruses"] = 0
        return empty

    by_genome = Counter()
    bp_by_genome = Counter()
    topology = Counter()
    taxonomy = Counter()
    scores = []
    score_labels = []
    lengths = []
    length_labels = []
    length_genomes = []
    with_amr = 0
    with_conj = 0
    with_tax = 0

    for row in rows:
        genome = _pretty(row.get("Genome")) or "unknown"
        contig = _pretty(row.get("Contig")) or "contig"
        topo = _pretty(row.get("Topology")) or "unknown"
        tax = _pretty(row.get("Taxonomy")) or "unclassified"
        score = _as_float(row.get("Score"))
        length = _as_int(row.get("Length")) or 0
        amr = _pretty(row.get("AMR genes"))
        conj = _pretty(row.get("Conjugation genes"))

        by_genome[genome] += 1
        bp_by_genome[genome] += max(length, 0)
        topology[topo] += 1
        taxonomy[tax] += 1
        if score is not None:
            scores.append(score)
            score_labels.append(f"{genome}:{contig}")
        lengths.append(max(length, 0))
        length_labels.append(contig)
        length_genomes.append(genome)
        if amr:
            with_amr += 1
        if conj:
            with_conj += 1
        if tax and tax.lower() not in ("unclassified", "na", "none"):
            with_tax += 1

    genomes = sorted(by_genome.keys())
    mean_score = round(sum(scores) / len(scores), 3) if scores else None
    ranked = sorted(
        zip(length_labels, lengths, length_genomes),
        key=lambda x: x[1],
        reverse=True,
    )[:40]

    # Taxonomy pie: top 8 + Other
    tax_items = taxonomy.most_common()
    if len(tax_items) > 8:
        top = tax_items[:7]
        other = sum(c for _, c in tax_items[7:])
        tax_labels = [k for k, _ in top] + ["Other"]
        tax_counts = [c for _, c in top] + [other]
    else:
        tax_labels = [k for k, _ in tax_items]
        tax_counts = [c for _, c in tax_items]

    summary = {
        "n_items": len(rows),
        "n_genomes": len(genomes),
        "mean_score": mean_score,
        "total_bp": int(sum(lengths)),
        "with_amr": with_amr,
        "with_conjugation": with_conj,
        "with_taxonomy": with_tax,
    }
    if kind == "plasmid":
        summary["n_plasmids"] = len(rows)
    else:
        summary["n_viruses"] = len(rows)

    return {
        "summary": summary,
        "by_genome": {
            "genomes": genomes,
            "counts": [by_genome[g] for g in genomes],
            "total_bp": [bp_by_genome[g] for g in genomes],
        },
        "scores": {"values": scores, "labels": score_labels},
        "lengths": {
            "contigs": [x[0] for x in ranked],
            "lengths": [x[1] for x in ranked],
            "genomes": [x[2] for x in ranked],
        },
        "topology": {
            "labels": list(topology.keys()),
            "counts": [topology[k] for k in topology.keys()],
        },
        "taxonomy": {"labels": tax_labels, "counts": tax_counts},
        "kind": kind,
    }


def plasmid_chart_payload(rows):
    return _chart_payload(rows, "plasmid")


def virus_chart_payload(rows):
    return _chart_payload(rows, "virus")
