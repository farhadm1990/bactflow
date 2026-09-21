#!/usr/bin/env python3
"""Parse CheckM QA tables and prune taxon Newick trees for the post-assembly UI."""

from __future__ import annotations

import os
import re
from collections import Counter

CHECKM_TABLE_COLUMNS = [
    "Bin Id",
    "Marker lineage",
    "# genomes",
    "# markers",
    "# marker sets",
    "0",
    "1",
    "2",
    "3",
    "4",
    "5+",
    "Completeness",
    "Contamination",
    "Strain heterogeneity",
]

_NUMERIC_TRAIL = 12
_RANK_RE = re.compile(r"^([a-z])__(.+)$")
_EXT_RE = re.compile(r"\.(fasta|fa|fna|faa|fsa|ffn)$", re.I)


def normalize_genome_id(value):
    text = os.path.basename(str(value or "").strip())
    text = _EXT_RE.sub("", text)
    return text


def _to_number(value):
    try:
        if "." in str(value):
            return float(value)
        return int(value)
    except (TypeError, ValueError):
        return value


def parse_checkm_lineage_table(text):
    """Parse CheckM qa pretty-table or TSV into a list of row dicts."""
    rows = []
    if not text:
        return rows
    lines = [ln.rstrip("\n") for ln in text.splitlines()]
    header_idx = None
    tab_header = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith("-") or stripped.startswith("["):
            continue
        if "Bin Id" in stripped and "Completeness" in stripped:
            header_idx = i
            if "\t" in line:
                tab_header = [c.strip() for c in line.split("\t")]
            break
    if header_idx is None:
        return rows

    if tab_header:
        for line in lines[header_idx + 1 :]:
            if not line.strip() or set(line.strip()) <= {"-", " "}:
                continue
            if line.strip().startswith("["):
                continue
            parts = [p.strip() for p in line.split("\t")]
            if len(parts) < 2:
                continue
            row = {}
            for key, val in zip(tab_header, parts):
                row[key] = _to_number(val) if key != "Bin Id" and key != "Marker lineage" else val
            if row.get("Bin Id"):
                rows.append(_normalize_row(row))
        return rows

    for line in lines[header_idx + 1 :]:
        stripped = line.strip()
        if not stripped or set(stripped) <= {"-", " "}:
            continue
        if stripped.startswith("[") or stripped.startswith("INFO"):
            continue
        parts = stripped.split()
        if len(parts) < _NUMERIC_TRAIL + 1:
            continue
        nums = parts[-_NUMERIC_TRAIL:]
        rest = parts[:-_NUMERIC_TRAIL]
        if not rest:
            continue
        try:
            [float(x) for x in nums]
        except ValueError:
            continue
        bin_id = rest[0]
        lineage = " ".join(rest[1:]) if len(rest) > 1 else ""
        values = [bin_id, lineage] + [_to_number(x) for x in nums]
        rows.append(dict(zip(CHECKM_TABLE_COLUMNS, values)))
    return rows


def _normalize_row(row):
    out = {}
    for col in CHECKM_TABLE_COLUMNS:
        out[col] = row.get(col, "")
    if not out["Bin Id"]:
        out["Bin Id"] = row.get("Bin Id") or row.get("Bin_Id") or ""
    return out


def classification_to_species(classification):
    ranks = [p.strip() for p in str(classification or "").split(";") if p.strip()]
    for prefix, suffix in (("s__", ""), ("g__", " sp."), ("f__", " sp."), ("o__", " sp.")):
        hit = next((p[3:] for p in reversed(ranks) if p.startswith(prefix) and p[3:]), "")
        if hit:
            return hit.replace("_", " ") + suffix
    return ""


def load_gtdb_species_map(out_dir):
    mapping = {}
    paths = []
    try:
        from gtdb_report import find_gtdb_summaries, salvage_gtdb_summaries
        salvage_gtdb_summaries(out_dir)
        paths = find_gtdb_summaries(out_dir)
    except Exception:
        paths = []
    if not paths:
        classify_dir = os.path.join(out_dir or "", "gtdbtk_out", "classify")
        if os.path.isdir(classify_dir):
            for name in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"):
                path = os.path.join(classify_dir, name)
                if os.path.isfile(path) and os.path.getsize(path) > 0:
                    paths.append(path)
    for path in paths:
        try:
            with open(path, encoding="utf-8", errors="replace") as handle:
                header = handle.readline()
                cols = [c.strip() for c in header.split("\t")]
                try:
                    g_idx = cols.index("user_genome")
                    c_idx = cols.index("classification")
                except ValueError:
                    continue
                for line in handle:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) <= max(g_idx, c_idx):
                        continue
                    genome = normalize_genome_id(parts[g_idx])
                    species = classification_to_species(parts[c_idx])
                    if genome and species:
                        mapping[genome] = species
        except OSError:
            continue
    return mapping


def lineage_to_species(lineage):
    text = str(lineage or "")
    match = re.search(r"([a-z])__([A-Za-z0-9_]+)", text)
    if match and match.group(2):
        return match.group(2).replace("_", " ")
    return ""


def species_for_bin(bin_id, lineage, gtdb_map):
    key = normalize_genome_id(bin_id)
    if key in gtdb_map:
        return gtdb_map[key]
    for cand, species in gtdb_map.items():
        if cand == key or cand.startswith(key) or key.startswith(cand):
            return species
    fallback = lineage_to_species(lineage)
    return fallback or key or str(bin_id)


def tip_key(label):
    text = str(label or "").strip().strip("'\"")
    if "|" in text:
        text = text.split("|", 1)[0]
    return normalize_genome_id(text)


def taxon_rank_label(raw):
    """Most specific CheckM/NCBI rank, e.g. o__Actinomycetales, or empty."""
    if not raw:
        return ""
    best = ""
    for chunk in re.split(r"[|;]", str(raw)):
        chunk = chunk.strip().strip("'\"")
        match = _RANK_RE.match(chunk)
        if not match:
            continue
        rest = match.group(2).strip()
        if not rest or rest.lower() in ("unresolved", "root"):
            continue
        best = f"{match.group(1)}__{rest}"
    return best


def species_from_tree_label(label):
    text = str(label or "")
    if "|" in text:
        tax = text.split("|", 1)[1]
        species = classification_to_species(tax)
        if species:
            return species
    ranked = taxon_rank_label(text)
    if ranked.startswith("s__"):
        return ranked[3:].replace("_", " ")
    return ""


def parse_newick(text):
    """Iterative Newick parser. Returns a node dict {name, length, children}."""
    s = re.sub(r"\[[^\]]*\]", "", text or "").strip()
    if not s or s[0] not in "(ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'_\"":
        raise ValueError("not a Newick tree")
    if s.endswith(";"):
        s = s[:-1]
    root = {"name": "", "length": None, "children": []}
    current = root
    stack = []
    buf = []
    in_quote = False
    quote_char = ""
    i = 0
    n = len(s)

    def flush(node):
        chunk = "".join(buf).strip()
        buf.clear()
        if not chunk:
            return
        name = chunk
        length = None
        if name[:1] in ("'", '"'):
            q = name[0]
            j = 1
            while j < len(name):
                if name[j] == q:
                    if j + 1 < len(name) and name[j + 1] == q:
                        j += 2
                        continue
                    inner = name[1:j].replace(q + q, q)
                    rest = name[j + 1 :]
                    name = inner
                    if rest.startswith(":"):
                        try:
                            length = float(rest[1:])
                        except ValueError:
                            length = None
                    break
                j += 1
            else:
                name = name.strip("'\"")
        else:
            if ":" in name:
                name, _, dist = name.rpartition(":")
                try:
                    length = float(dist)
                except ValueError:
                    length = None
        node["name"] = name
        node["length"] = length

    while i < n:
        c = s[i]
        if in_quote:
            buf.append(c)
            if c == quote_char:
                in_quote = False
            i += 1
            continue
        if c in ("'", '"'):
            in_quote = True
            quote_char = c
            buf.append(c)
            i += 1
            continue
        if c == "(":
            child = {"name": "", "length": None, "children": []}
            current["children"].append(child)
            stack.append(current)
            current = child
            i += 1
        elif c == ",":
            flush(current)
            parent = stack[-1]
            child = {"name": "", "length": None, "children": []}
            parent["children"].append(child)
            current = child
            i += 1
        elif c == ")":
            flush(current)
            current = stack.pop()
            i += 1
        else:
            buf.append(c)
            i += 1
    flush(current)
    return root


def iter_nodes(root):
    stack = [root]
    while stack:
        node = stack.pop()
        yield node
        stack.extend(node["children"])


def count_tips(root):
    return sum(1 for node in iter_nodes(root) if not node["children"])


def _ids_index(bin_ids):
    return {normalize_genome_id(x) for x in bin_ids if x}


def tip_matches(label, id_set):
    key = tip_key(label)
    if not key:
        return False
    if key in id_set:
        return True
    for item in id_set:
        if key.startswith(item) or item.startswith(key):
            return True
    return False


def prune_to_ids(root, bin_ids):
    id_set = _ids_index(bin_ids)
    if not id_set or root is None:
        return None
    nodes = list(iter_nodes(root))
    parent = {id(root): None}
    for node in nodes:
        for child in node["children"]:
            parent[id(child)] = node
    keep = set()
    for node in nodes:
        if node["children"]:
            continue
        if not tip_matches(node.get("name") or "", id_set):
            continue
        cur = node
        while cur is not None and id(cur) not in keep:
            keep.add(id(cur))
            cur = parent.get(id(cur))
    if id(root) not in keep:
        return None
    for node in nodes:
        node["children"] = [c for c in node["children"] if id(c) in keep]
    return root


def collapse_unary(node, is_root=True):
    node["children"] = [collapse_unary(child, False) for child in node["children"]]
    named = taxon_rank_label(node.get("name") or "")
    if not is_root and len(node["children"]) == 1 and not named:
        child = node["children"][0]
        child["length"] = (child.get("length") or 0.0) + (node.get("length") or 0.0)
        return child
    if is_root and len(node["children"]) == 1 and not named:
        return node["children"][0]
    return node


def quote_newick_name(name):
    text = str(name or "")
    if not text:
        return ""
    if re.search(r"[\s(),:;'\[\]]", text):
        return "'" + text.replace("'", "''") + "'"
    return text


def format_branch_length(value):
    if value is None:
        return ""
    text = f"{float(value):.8f}".rstrip("0").rstrip(".")
    if not text or text == "-0":
        text = "0"
    return f":{text}"


def to_newick(node):
    def rec(n):
        label = quote_newick_name(n.get("name") or "")
        length = format_branch_length(n.get("length"))
        children = n.get("children") or []
        if children:
            return "(" + ",".join(rec(c) for c in children) + ")" + label + length
        return label + length

    return rec(node) + ";"


def relabel_tree(node, species_map, bin_ids):
    used = []

    def lookup_species(key, raw_name):
        if key in species_map:
            return species_map[key]
        for cand, species in species_map.items():
            if cand == key or (cand and key and (cand.startswith(key) or key.startswith(cand))):
                return species
        return species_from_tree_label(raw_name) or key

    def walk(n):
        if n.get("children"):
            n["name"] = taxon_rank_label(n.get("name") or "")
            for child in n["children"]:
                walk(child)
            return
        key = tip_key(n.get("name") or "")
        n["name"] = lookup_species(key, n.get("name") or "") or key
        used.append((n, key))

    walk(node)
    counts = Counter(n["name"] for n, _ in used)
    for n, key in used:
        if counts[n["name"]] > 1 and key and key not in n["name"]:
            n["name"] = f"{n['name']} ({key})"
    return node


def build_species_map(rows, gtdb_map):
    mapping = {}
    for row in rows:
        bin_id = normalize_genome_id(row.get("Bin Id"))
        if not bin_id:
            continue
        mapping[bin_id] = species_for_bin(bin_id, row.get("Marker lineage"), gtdb_map)
    return mapping


def prepare_checkm_tree(newick_text, rows, gtdb_map):
    text = (newick_text or "").strip()
    if not text or not text.startswith("(") and not text[:1].isalnum():
        return None
    if not text.startswith("(") and "Bin Id" in text:
        return None
    try:
        tree = parse_newick(text)
    except (ValueError, IndexError):
        return None
    bin_ids = [r.get("Bin Id") for r in rows if r.get("Bin Id")]
    species_map = build_species_map(rows, gtdb_map)
    n_tips = count_tips(tree)
    pruned = prune_to_ids(tree, bin_ids) if bin_ids else None
    if pruned is None:
        if n_tips > 40:
            return None
        pruned = tree
    pruned = collapse_unary(pruned)
    if count_tips(pruned) < 1:
        return None
    relabel_tree(pruned, species_map, bin_ids)
    return to_newick(pruned)


def read_text_if_nonempty(path, max_bytes=80_000_000):
    if not path or not os.path.isfile(path):
        return ""
    try:
        size = os.path.getsize(path)
        if size <= 0 or size > max_bytes:
            return ""
        with open(path, encoding="utf-8", errors="replace") as handle:
            return handle.read().strip()
    except OSError:
        return ""


def find_checkm_files(out_dir):
    if not out_dir:
        return "", "", ""
    candidates = [
        os.path.join(out_dir, "checkm_out"),
        os.path.join(out_dir, "checkm_lineage"),
        out_dir,
    ]
    lineage = ""
    taxon = ""
    genome = ""
    for folder in candidates:
        lin = os.path.join(folder, "checkm_lineage.txt")
        if not os.path.isfile(lin) or os.path.getsize(lin) <= 0:
            continue
        lineage = lin
        tax = os.path.join(folder, "taxon_tree.newick")
        gen = os.path.join(folder, "genome_tree.tree")
        gen2 = os.path.join(folder, "genome_tree.newick")
        if os.path.isfile(tax) and os.path.getsize(tax) > 0:
            taxon = tax
        if os.path.isfile(gen) and os.path.getsize(gen) > 0:
            genome = gen
        elif os.path.isfile(gen2) and os.path.getsize(gen2) > 0:
            genome = gen2
        break
    return lineage, taxon, genome
