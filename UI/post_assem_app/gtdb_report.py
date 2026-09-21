#!/usr/bin/env python3
"""Find GTDB-Tk summary tables and build taxonomy + abundance views."""

from __future__ import annotations

import os
import glob
import shutil

import pandas as pd

RANK_PREFIXES = [
    ("Domain", "d__"),
    ("Phylum", "p__"),
    ("Class", "c__"),
    ("Order", "o__"),
    ("Family", "f__"),
    ("Genus", "g__"),
    ("Species", "s__"),
]

ANI_COLS = (
    "closest_genome_ani",
    "fastani_ani",
    "ani",
    "closest_ani",
    "closest_genome_af",
)

SUMMARY_NAMES = (
    "gtdbtk.bac120.summary.tsv",
    "gtdbtk.ar53.summary.tsv",
)


def _pretty_taxon(raw):
    text = str(raw or "").strip()
    if not text or text in (".", "nan", "None"):
        return "unclassified"
    if len(text) > 3 and text[1:3] == "__":
        rest = text[3:].replace("_", " ").strip()
        return rest or "unclassified"
    return text.replace("_", " ")


def extract_rank(classification, prefix):
    for part in str(classification or "").split(";"):
        part = part.strip()
        if part.startswith(prefix):
            rest = part[len(prefix) :].strip()
            if rest:
                return part
    return ""


def _is_summary_file(path):
    if not os.path.isfile(path) or os.path.getsize(path) <= 0:
        return False
    name = os.path.basename(path).lower()
    if "markers" in name:
        return False
    return name.endswith(".summary.tsv")


def _collect_from_dir(root, recursive, seen, found):
    if not root or not os.path.isdir(root):
        return
    if not recursive:
        try:
            names = os.listdir(root)
        except OSError:
            return
        for name in names:
            path = os.path.join(root, name)
            if not _is_summary_file(path):
                continue
            key = os.path.realpath(path)
            if key in seen:
                continue
            seen.add(key)
            found.append(path)
        return
    try:
        walker = os.walk(root)
    except OSError:
        return
    for dirpath, _dirnames, names in walker:
        for name in names:
            path = os.path.join(dirpath, name)
            if not _is_summary_file(path):
                continue
            key = os.path.realpath(path)
            if key in seen:
                continue
            seen.add(key)
            found.append(path)


def find_gtdb_summaries(out_dir):
    if not out_dir:
        return []
    found = []
    seen = set()
    _collect_from_dir(os.path.join(out_dir, "gtdbtk_out", "classify"), False, seen, found)
    _collect_from_dir(os.path.join(out_dir, "gtdbtk_out"), True, seen, found)
    _collect_from_dir(out_dir, False, seen, found)
    return found


def _copy_into_classify(out_dir, sources):
    dest_dir = os.path.join(out_dir, "gtdbtk_out", "classify")
    copied = []
    try:
        os.makedirs(dest_dir, exist_ok=True)
    except OSError:
        return copied
    for path in sources:
        if not _is_summary_file(path):
            continue
        target = os.path.join(dest_dir, os.path.basename(path))
        try:
            if os.path.realpath(path) != os.path.realpath(target):
                shutil.copy2(path, target)
            copied.append(target)
        except OSError:
            continue
    return copied


def _iter_nextflow_task_dirs(work_root):
    if not work_root or not os.path.isdir(work_root):
        return
    try:
        level1 = os.listdir(work_root)
    except OSError:
        return
    for hash1 in level1:
        p1 = os.path.join(work_root, hash1)
        if not os.path.isdir(p1):
            continue
        try:
            level2 = os.listdir(p1)
        except OSError:
            continue
        for hash2 in level2:
            p2 = os.path.join(p1, hash2)
            if os.path.isdir(p2):
                yield p2


def salvage_gtdb_summaries(out_dir, extra_roots=None):
    """Copy GTDB-Tk summaries into out_dir/gtdbtk_out/classify, including from Nextflow work."""
    if not out_dir:
        return []
    existing = find_gtdb_summaries(out_dir)
    if existing:
        _copy_into_classify(out_dir, existing)
        return find_gtdb_summaries(out_dir)

    roots = []
    if extra_roots:
        roots.extend(extra_roots)
    roots.extend([
        os.environ.get("NXF_WORK") or "",
        os.path.join(out_dir, ".nextflow-work"),
        os.path.join(out_dir, "work"),
        os.path.join(os.path.dirname(os.path.abspath(out_dir)), "work"),
        os.path.join(os.getcwd(), "work"),
        os.path.join(os.getcwd(), "bactflow_out", ".nextflow-work"),
        os.path.join(os.environ.get("HOME") or "", "work"),
    ])

    found = []
    seen = set()
    for root in roots:
        if not root or not os.path.isdir(root):
            continue
        patterns = [
            os.path.join(root, "*", "*", "*.summary.tsv"),
            os.path.join(root, "*", "*", "gtdbtk_out", "classify", "*.summary.tsv"),
            os.path.join(root, "*", "*", "classify", "*.summary.tsv"),
        ]
        for pat in patterns:
            for cand in glob.glob(pat):
                if not _is_summary_file(cand):
                    continue
                key = os.path.realpath(cand)
                if key in seen:
                    continue
                seen.add(key)
                found.append(cand)
        if found:
            break
    if found:
        _copy_into_classify(out_dir, found)
    return find_gtdb_summaries(out_dir)


def load_gtdb_taxonomy(out_dir, extra_roots=None):
    salvage_gtdb_summaries(out_dir, extra_roots=extra_roots)
    frames = []
    for path in find_gtdb_summaries(out_dir):
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception:
            continue
        if df.empty or "user_genome" not in df.columns or "classification" not in df.columns:
            continue
        ani_col = next((c for c in ANI_COLS if c in df.columns), None)
        keep = ["user_genome", "classification"]
        if ani_col:
            keep.append(ani_col)
        part = df[keep].copy()
        if ani_col and ani_col != "closest_genome_ani":
            part = part.rename(columns={ani_col: "closest_genome_ani"})
        elif "closest_genome_ani" not in part.columns:
            part["closest_genome_ani"] = ""
        frames.append(part)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    out = out.drop_duplicates(subset=["user_genome"], keep="first")
    return out


def taxonomy_display_rows(df):
    rows = []
    for rec in df.to_dict(orient="records"):
        classif = rec.get("classification")
        try:
            if classif is None or pd.isna(classif):
                classif = ""
            else:
                classif = str(classif)
        except (TypeError, ValueError):
            classif = str(classif or "")
        ani = rec.get("closest_genome_ani")
        try:
            if ani is None or pd.isna(ani):
                ani = ""
        except (TypeError, ValueError):
            if ani is None:
                ani = ""
        rows.append({
            "Genome": rec.get("user_genome") or "",
            "Domain": _pretty_taxon(extract_rank(classif, "d__")),
            "Phylum": _pretty_taxon(extract_rank(classif, "p__")),
            "Class": _pretty_taxon(extract_rank(classif, "c__")),
            "Order": _pretty_taxon(extract_rank(classif, "o__")),
            "Family": _pretty_taxon(extract_rank(classif, "f__")),
            "Genus": _pretty_taxon(extract_rank(classif, "g__")),
            "Species": _pretty_taxon(extract_rank(classif, "s__")),
            "Classification": classif,
            "Closest ANI": ani,
        })
    return rows


def abundance_display_rows(df):
    total = max(len(df), 1)
    rows = []
    for rank_name, prefix in RANK_PREFIXES:
        counts = {}
        for classif in df.get("classification", []):
            key = extract_rank(classif, prefix) or f"{prefix}unclassified"
            counts[key] = counts.get(key, 0) + 1
        for taxon, count in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
            rows.append({
                "Rank": rank_name,
                "Taxon": _pretty_taxon(taxon),
                "Count": int(count),
                "Percent": round(100.0 * count / total, 1),
            })
    return rows


def find_gtdb_trees(out_dir):
    if not out_dir:
        return []
    found = []
    seen = set()
    roots = [
        os.path.join(out_dir, "gtdbtk_out", "classify"),
        os.path.join(out_dir, "gtdbtk_out"),
    ]
    patterns = ("*.classify.tree", "*.decorated.tree")
    for root in roots:
        if not os.path.isdir(root):
            continue
        for pat in patterns:
            for path in glob.glob(os.path.join(root, pat)):
                if not os.path.isfile(path) or os.path.getsize(path) <= 0:
                    continue
                key = os.path.realpath(path)
                if key in seen:
                    continue
                seen.add(key)
                found.append(path)
        if found:
            break
    found.sort(key=lambda p: (0 if "classify.tree" in os.path.basename(p) else 1, os.path.basename(p)))
    return found


def prepare_gtdb_tree(out_dir, extra_roots=None):
    """Prune GTDB-Tk classify trees to user genomes and relabel tips with species."""
    df = load_gtdb_taxonomy(out_dir, extra_roots=extra_roots)
    if df is None or df.empty:
        return "", ""
    try:
        from checkm_report import (
            classification_to_species,
            collapse_unary,
            count_tips,
            normalize_genome_id,
            parse_newick,
            prune_to_ids,
            read_text_if_nonempty,
            relabel_tree,
            to_newick,
        )
    except Exception:
        return "", ""
    genome_ids = []
    species_map = {}
    for rec in df.to_dict(orient="records"):
        genome = rec.get("user_genome") or ""
        if not genome:
            continue
        genome_ids.append(genome)
        key = normalize_genome_id(genome)
        species = classification_to_species(rec.get("classification") or "")
        if key and species:
            species_map[key] = species
    if not genome_ids:
        return "", ""
    for path in find_gtdb_trees(out_dir):
        raw = read_text_if_nonempty(path)
        if not raw or not raw.lstrip().startswith("("):
            continue
        try:
            tree = parse_newick(raw)
        except (ValueError, IndexError):
            continue
        pruned = prune_to_ids(tree, genome_ids)
        if pruned is None:
            continue
        pruned = collapse_unary(pruned)
        if count_tips(pruned) < 1:
            continue
        relabel_tree(pruned, species_map, genome_ids)
        return to_newick(pruned), os.path.basename(path)
    return "", ""
