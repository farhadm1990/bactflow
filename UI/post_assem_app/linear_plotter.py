"""Interactive linear genome maps from Bakta GFF / GenBank.

The layout follows DNA Features Viewer (Zulkower & Rosser, Bioinformatics 2020):
annotated genes sit on the forward strand above a backbone and on the reverse
strand below, colored by feature type. The figure is Plotly so the whole
genome stays visible with hover, zoom, and pan.
"""

from __future__ import annotations

import glob
import os

from Bio import SeqIO

GBK_EXTS = {".gbk", ".gbff", ".gb"}
GFF_EXTS = {".gff", ".gff3"}
ANNOT_EXTS = GBK_EXTS | GFF_EXTS
SKIP_TYPES = {"source", "region", "comment"}

TYPE_COLORS = {
    "cds": "#3B82F6",
    "trna": "#F59E0B",
    "rrna": "#EF4444",
    "tmrna": "#14B8A6",
    "ncrna": "#22C55E",
    "repeat_region": "#A855F7",
    "gene": "#64748B",
    "orig": "#EAB308",
    "oric": "#EAB308",
    "gap": "#94A3B8",
}
DEFAULT_COLOR = "#64748B"


def _first(attrs, *keys):
    for key in keys:
        val = attrs.get(key)
        if isinstance(val, list):
            val = val[0] if val else ""
        if val:
            return str(val)
    return ""


def _parse_gff_attrs(col):
    attrs = {}
    for part in (col or "").split(";"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        attrs[key.strip()] = value.strip()
    return attrs


def _companion_fasta(annot_path):
    stem, _ = os.path.splitext(annot_path)
    for ext in (".fna", ".fa", ".fasta"):
        candidate = stem + ext
        if os.path.isfile(candidate):
            return candidate
    parent = os.path.dirname(annot_path)
    for ext in (".fna", ".fa", ".fasta"):
        hits = glob.glob(os.path.join(parent, f"*{ext}"))
        if hits:
            return hits[0]
    return None


def _ensure_contig(store, order, cid, length=0):
    if cid not in store:
        order.append(cid)
        store[cid] = {"id": cid, "length": max(int(length), 0), "features": []}
    elif length:
        store[cid]["length"] = max(store[cid]["length"], int(length))
    return store[cid]


def _add_feature(contig, ftype, start, end, strand, name, product, locus):
    if ftype.lower() in SKIP_TYPES:
        return
    start = max(int(start), 1)
    end = max(int(end), start)
    contig["length"] = max(contig["length"], end)
    contig["features"].append({
        "type": ftype,
        "start": start,
        "end": end,
        "strand": "-" if str(strand) in {"-", "-1", "-1.0"} else "+",
        "name": name or ftype,
        "product": product or "",
        "locus": locus or "",
    })


def parse_gff(path):
    store, order = {}, []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                break
            if line.startswith("##sequence-region"):
                parts = line.split()
                if len(parts) >= 4:
                    _ensure_contig(store, order, parts[1], int(parts[3]) - int(parts[2]) + 1)
                continue
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 7:
                continue
            attrs = _parse_gff_attrs(cols[8] if len(cols) > 8 else "")
            contig = _ensure_contig(store, order, cols[0], int(cols[4]))
            _add_feature(
                contig,
                cols[2],
                cols[3],
                cols[4],
                cols[6],
                _first(attrs, "Name", "gene", "locus_tag", "ID"),
                _first(attrs, "product", "Product", "Note"),
                _first(attrs, "locus_tag", "ID"),
            )

    fasta = _companion_fasta(path)
    if fasta:
        for rec in SeqIO.parse(fasta, "fasta"):
            cid = rec.id.split()[0]
            _ensure_contig(store, order, cid, len(rec.seq))
    return [store[cid] for cid in order if cid in store]


def parse_gbk(path):
    store, order = {}, []
    for rec in SeqIO.parse(path, "genbank"):
        cid = rec.id.split()[0] if rec.id else rec.name or os.path.basename(path)
        contig = _ensure_contig(store, order, cid, len(rec.seq))
        for feat in rec.features:
            if feat.location is None:
                continue
            start = int(feat.location.start) + 1
            end = int(feat.location.end)
            strand = "+" if feat.location.strand != -1 else "-"
            quals = feat.qualifiers or {}
            _add_feature(
                contig,
                feat.type,
                start,
                end,
                strand,
                _first(quals, "gene", "locus_tag", "label"),
                _first(quals, "product", "note"),
                _first(quals, "locus_tag"),
            )
    return [store[cid] for cid in order if cid in store]


def collect_annotation_files(root):
    chosen = []
    for dirpath, _dirnames, filenames in os.walk(root):
        by_ext = {ext: [] for ext in ANNOT_EXTS}
        for name in filenames:
            ext = os.path.splitext(name)[1].lower()
            if ext in ANNOT_EXTS:
                by_ext[ext].append(os.path.join(dirpath, name))
        picked = None
        for ext in (".gbff", ".gbk", ".gb", ".gff3", ".gff"):
            if by_ext.get(ext):
                picked = sorted(by_ext[ext])[0]
                break
        if picked:
            chosen.append(picked)
    return chosen


def load_genome(path):
    ext = os.path.splitext(path)[1].lower()
    name = os.path.splitext(os.path.basename(path))[0]
    if ext in GFF_EXTS:
        contigs = parse_gff(path)
    elif ext in GBK_EXTS:
        contigs = parse_gbk(path)
    else:
        raise ValueError(f"Unsupported annotation file: {path}")
    return {"name": name, "path": path, "contigs": contigs}


def _wanted_types(feature_types):
    raw = feature_types if isinstance(feature_types, (list, tuple)) else str(feature_types or "cds").split(",")
    wanted = {str(t).strip().lower() for t in raw if str(t).strip()}
    return wanted or {"cds"}


def _pack_lanes(features):
    lanes = []
    for feat in features:
        placed = False
        for idx, last_end in enumerate(lanes):
            if feat["start"] > last_end + 8:
                lanes[idx] = feat["end"]
                feat["lane"] = idx
                placed = True
                break
        if not placed:
            feat["lane"] = len(lanes)
            lanes.append(feat["end"])
    return max(len(lanes), 1)


def _type_color(ftype, reverse=False):
    color = TYPE_COLORS.get(ftype.lower(), DEFAULT_COLOR)
    if not reverse:
        return color
    raw = color.lstrip("#")
    if len(raw) != 6:
        return color
    r, g, b = (int(raw[i:i + 2], 16) for i in (0, 2, 4))
    return f"#{max(0, int(r * 0.72)):02x}{max(0, int(g * 0.72)):02x}{max(0, int(b * 0.72)):02x}"


def _nice_tick(length):
    length = max(int(length), 1)
    target = max(length / 12.0, 1.0)
    for step in (500, 1000, 2000, 5000, 10000, 20000, 25000, 50000, 100000, 200000, 250000, 500000, 1000000, 2000000, 5000000):
        if step >= target:
            return step
    return 10000000


def _fmt_bp(value):
    value = int(value)
    if abs(value) >= 1_000_000:
        return f"{value / 1e6:.2f} Mb"
    if abs(value) >= 10_000:
        return f"{value / 1e3:.1f} kb"
    return f"{value} bp"


def build_linear_figure(bakta_dir, feature_types=None, genome="", contig=""):
    files = collect_annotation_files(bakta_dir)
    if not files:
        return None, {"error": "No .gbk/.gbff/.gff files found in Bakta output.", "genomes": [], "contigs": []}

    genomes = []
    for path in files:
        try:
            genomes.append(load_genome(path))
        except Exception:
            continue
    if not genomes:
        return None, {"error": "Could not parse Bakta annotation files.", "genomes": [], "contigs": []}

    names = [g["name"] for g in genomes]
    chosen = next((g for g in genomes if g["name"] == genome), genomes[0])
    contig_ids = [c["id"] for c in chosen["contigs"]]
    wanted = _wanted_types(feature_types)

    if contig and contig != "__all__":
        pieces = [c for c in chosen["contigs"] if c["id"] == contig]
        if not pieces:
            pieces = list(chosen["contigs"])
            contig = "__all__"
    else:
        pieces = list(chosen["contigs"])
        contig = contig or "__all__"

    if not pieces:
        return None, {
            "error": "No contigs in the selected genome.",
            "genomes": names,
            "contigs": contig_ids,
            "selected": {"genome": chosen["name"], "contig": contig},
        }

    rows = []
    offset = 0
    boundaries = []
    total = 0
    for piece in pieces:
        length = max(int(piece.get("length") or 0), 1)
        boundaries.append({"id": piece["id"], "start": offset, "end": offset + length, "length": length})
        for feat in piece.get("features") or []:
            if feat["type"].lower() not in wanted:
                continue
            rows.append({
                **feat,
                "contig": piece["id"],
                "abs_start": offset + feat["start"],
                "abs_end": offset + feat["end"],
            })
        offset += length
        total = offset

    fwd = [f for f in rows if f["strand"] != "-"]
    rev = [f for f in rows if f["strand"] == "-"]
    fwd_lanes = _pack_lanes(fwd)
    rev_lanes = _pack_lanes(rev)
    lane_h = 0.42
    y_max = 0.55 + fwd_lanes * lane_h
    y_min = -(0.55 + rev_lanes * lane_h)

    features = []
    seen_types = []
    for feat in fwd + rev:
        ftype = feat["type"].lower()
        if ftype not in seen_types:
            seen_types.append(ftype)
        features.append({
            "s": int(feat["abs_start"]),
            "e": int(feat["abs_end"]),
            "lane": int(feat.get("lane") or 0),
            "strand": feat["strand"],
            "t": ftype,
            "n": (feat["name"] or ftype)[:80],
            "p": (feat["product"] or "")[:120],
            "c": feat["contig"],
        })

    title_contig = "all contigs" if contig == "__all__" or len(pieces) > 1 else pieces[0]["id"]
    return {"canvas": True}, {
        "genomes": names,
        "contigs": contig_ids,
        "selected": {"genome": chosen["name"], "contig": contig},
        "feature_count": len(rows),
        "length": total,
        "length_label": _fmt_bp(total),
        "title": f"{chosen['name']} · {title_contig}",
        "fwd_lanes": fwd_lanes,
        "rev_lanes": rev_lanes,
        "types": seen_types,
        "features": features,
        "boundaries": [{"id": b["id"], "start": b["start"]} for b in boundaries[1:]],
    }
