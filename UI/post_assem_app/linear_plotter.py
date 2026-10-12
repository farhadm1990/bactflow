"""Interactive linear genome maps from Bakta GFF / GenBank.

The layout follows DNA Features Viewer (Zulkower & Rosser, Bioinformatics 2020):
annotated genes sit on the forward strand above a backbone and on the reverse
strand below, colored by feature type. The figure is Plotly so the whole
genome stays visible with hover, zoom, and pan.
"""

from __future__ import annotations

import glob
import os
import re
from io import StringIO

from Bio import SeqIO

GBK_EXTS = {".gbk", ".gbff", ".gb"}
GFF_EXTS = {".gff", ".gff3"}
ANNOT_EXTS = GBK_EXTS | GFF_EXTS
SKIP_TYPES = {"source", "region", "comment"}

TYPE_COLORS = {
    "cds": "#E07A5F",
    "trna": "#E9C46A",
    "rrna": "#E76F51",
    "tmrna": "#2A9D8F",
    "ncrna": "#6A994E",
    "repeat_region": "#4C6EF5",
    "crispr": "#9B5DE5",
    "crispr-repeat": "#C77DFF",
    "crispr-spacer": "#B8C0FF",
    "ncrna-region": "#A7C957",
    "gene": "#6D6875",
    "orig": "#F4A261",
    "oric": "#F4A261",
    "gap": "#94A3B8",
}
DEFAULT_COLOR = "#6D6875"


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


def _norm_ftype(ftype, name="", product=""):
    raw = str(ftype or "")
    t = raw.lower().replace("_", "-")
    blob = f"{name or ''} {product or ''} {raw}".lower()
    if t in {"crispr-spacer", "crisprspacer"} or "crispr spacer" in blob:
        return "crispr-spacer"
    if t in {"crispr-repeat", "crisprrepeat"} or "crispr repeat" in blob:
        return "crispr-repeat"
    if t == "crispr" or (t == "repeat-region" and "crispr" in blob):
        if "spacer" in blob:
            return "crispr-spacer"
        if "repeat" in blob and "array" not in blob:
            return "crispr-repeat"
        return "crispr"
    if t == "repeat-region":
        return "repeat_region"
    return raw.lower()


def _add_feature(contig, ftype, start, end, strand, name, product, locus):
    if ftype.lower() in SKIP_TYPES:
        return
    start = max(int(start), 1)
    end = max(int(end), start)
    contig["length"] = max(contig["length"], end)
    contig["features"].append({
        "type": _norm_ftype(ftype, name, product),
        "start": start,
        "end": end,
        "strand": "-" if str(strand) in {"-", "-1", "-1.0"} else "+",
        "name": name or ftype,
        "product": product or "",
        "locus": locus or "",
    })


def parse_gff(path):
    store, order = {}, []
    fasta_chunks = []
    in_fasta = False
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("##FASTA"):
                in_fasta = True
                continue
            if in_fasta:
                fasta_chunks.append(line)
                continue
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

    def _attach_seqs(records):
        for rec in records:
            cid = rec.id.split()[0]
            contig = _ensure_contig(store, order, cid, len(rec.seq))
            if not contig.get("seq"):
                contig["seq"] = str(rec.seq).upper()

    if fasta_chunks:
        _attach_seqs(SeqIO.parse(StringIO("".join(fasta_chunks)), "fasta"))
    fasta = _companion_fasta(path)
    if fasta:
        _attach_seqs(SeqIO.parse(fasta, "fasta"))
    return [store[cid] for cid in order if cid in store]


def parse_gbk(path):
    store, order = {}, []
    for rec in SeqIO.parse(path, "genbank"):
        cid = rec.id.split()[0] if rec.id else rec.name or os.path.basename(path)
        contig = _ensure_contig(store, order, cid, len(rec.seq))
        contig["seq"] = str(rec.seq).upper()
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
    rank = {".gbff": 0, ".gbk": 1, ".gb": 2, ".gff3": 3, ".gff": 4}
    by_stem = {}
    for dirpath, _dirnames, filenames in os.walk(root):
        for name in filenames:
            ext = os.path.splitext(name)[1].lower()
            if ext not in ANNOT_EXTS:
                continue
            stem = os.path.splitext(name)[0]
            key = os.path.join(dirpath, stem)
            prev = by_stem.get(key)
            prev_rank = rank.get(os.path.splitext(prev)[1].lower(), 9) if prev else 9
            if prev is None or rank.get(ext, 9) < prev_rank:
                by_stem[key] = os.path.join(dirpath, name)
    return sorted(by_stem.values())


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


_GENOME_CACHE = {}


def load_genome_cached(path):
    path = os.path.abspath(path)
    try:
        stamp = (os.path.getmtime(path), os.path.getsize(path))
    except OSError:
        stamp = (0, 0)
    hit = _GENOME_CACHE.get(path)
    if hit and hit[0] == stamp:
        return hit[1]
    genome = load_genome(path)
    _GENOME_CACHE[path] = (stamp, genome)
    return genome


def _wanted_types(feature_types):
    raw = feature_types if isinstance(feature_types, (list, tuple)) else str(feature_types or "cds").split(",")
    wanted = {str(t).strip().lower() for t in raw if str(t).strip()}
    if not wanted:
        wanted = {"cds"}
    if wanted & {"ncrna", "ncrna-region"}:
        wanted.update({"ncrna", "ncrna-region"})
    if "crispr-repeat" in wanted:
        wanted.update({"crispr-repeat", "crispr-spacer"})
    if wanted & {"repeat", "repeat_region"}:
        wanted.add("repeat_region")
    return wanted


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


_HYPO_RE = re.compile(
    r"hypothetical|uncharacterized|unknown function|putative protein|unannotated",
    re.I,
)


def _is_hypothetical(feat):
    if str(feat.get("type") or "").lower() != "cds":
        return False
    product = str(feat.get("product") or "").strip()
    name = str(feat.get("name") or "").strip()
    if not product:
        return True
    if _HYPO_RE.search(product):
        return True
    return name.lower() in {"cds", "hypothetical protein"}


def _parse_genome_names(genome):
    if isinstance(genome, (list, tuple)):
        return [str(item).strip() for item in genome if str(item).strip()]
    return [part.strip() for part in str(genome or "").split(",") if part.strip()]


def _track_from_genome(chosen, contig, wanted):
    if contig and contig != "__all__":
        pieces = [c for c in chosen["contigs"] if c["id"] == contig]
        if not pieces:
            pieces = list(chosen["contigs"])
            contig = "__all__"
    else:
        pieces = list(chosen["contigs"])
        contig = contig or "__all__"

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
    features = []
    seen_types = []
    hypo_count = 0
    for feat in fwd + rev:
        ftype = feat["type"].lower()
        if ftype not in seen_types:
            seen_types.append(ftype)
        hypo = _is_hypothetical(feat)
        if hypo:
            hypo_count += 1
        features.append({
            "s": int(feat["abs_start"]),
            "e": int(feat["abs_end"]),
            "lane": int(feat.get("lane") or 0),
            "strand": feat["strand"],
            "t": ftype,
            "n": (feat["name"] or ftype)[:80],
            "p": (feat["product"] or "")[:120],
            "c": feat["contig"],
            "h": 1 if hypo else 0,
        })
    title_contig = "all contigs" if contig == "__all__" or len(pieces) > 1 else pieces[0]["id"]
    return {
        "name": chosen["name"],
        "contig": contig,
        "length": total,
        "length_label": _fmt_bp(total),
        "fwd_lanes": fwd_lanes,
        "rev_lanes": rev_lanes,
        "types": seen_types,
        "features": features,
        "boundaries": boundaries,
        "title": f"{chosen['name']} · {title_contig}",
        "hypo_count": hypo_count,
        "feature_count": len(features),
    }


def build_linear_figure(bakta_dir, feature_types=None, genome="", contig=""):
    files = collect_annotation_files(bakta_dir)
    if not files:
        return None, {"error": "No .gbk/.gbff/.gff files found in Bakta output.", "genomes": [], "contigs": []}

    genomes = []
    for path in files:
        try:
            genomes.append(load_genome_cached(path))
        except Exception:
            continue
    if not genomes:
        return None, {"error": "Could not parse Bakta annotation files.", "genomes": [], "contigs": []}

    names = [g["name"] for g in genomes]
    wanted = _wanted_types(feature_types)
    picked_names = _parse_genome_names(genome)
    if picked_names:
        chosen_list = [g for g in genomes if g["name"] in picked_names]
        order = {name: idx for idx, name in enumerate(picked_names)}
        chosen_list.sort(key=lambda item: order.get(item["name"], 999))
    else:
        chosen_list = [genomes[0]]
    if not chosen_list:
        chosen_list = [genomes[0]]

    tracks = []
    contig_ids = []
    for chosen in chosen_list:
        track = _track_from_genome(chosen, contig, wanted)
        tracks.append(track)
        for cid in [c["id"] for c in chosen["contigs"]]:
            if cid not in contig_ids:
                contig_ids.append(cid)

    if not tracks or not any(t["length"] for t in tracks):
        return None, {
            "error": "No contigs in the selected genome.",
            "genomes": names,
            "contigs": contig_ids,
            "selected": {"genomes": [t["name"] for t in tracks], "contig": contig or "__all__"},
        }

    types = []
    for track in tracks:
        for ftype in track["types"]:
            if ftype not in types:
                types.append(ftype)
    length = max(int(track["length"] or 1) for track in tracks)
    selected_names = [track["name"] for track in tracks]
    title = tracks[0]["title"] if len(tracks) == 1 else f"{len(tracks)} genomes"
    return {"canvas": True}, {
        "genomes": names,
        "contigs": contig_ids,
        "selected": {
            "genome": selected_names[0],
            "genomes": selected_names,
            "contig": tracks[0]["contig"],
        },
        "feature_count": sum(t["feature_count"] for t in tracks),
        "hypo_count": sum(t["hypo_count"] for t in tracks),
        "length": length,
        "length_label": _fmt_bp(length),
        "title": title,
        "fwd_lanes": tracks[0]["fwd_lanes"],
        "rev_lanes": tracks[0]["rev_lanes"],
        "types": types,
        "features": tracks[0]["features"],
        "boundaries": tracks[0]["boundaries"],
        "tracks": tracks,
    }


def _select_genome_pieces(bakta_dir, genome="", contig=""):
    files = collect_annotation_files(bakta_dir)
    if not files:
        return None, None, []
    if genome:
        named = [path for path in files if os.path.splitext(os.path.basename(path))[0] == genome]
        if named:
            files = named
    genomes = []
    for path in files:
        try:
            genomes.append(load_genome_cached(path))
        except Exception:
            continue
    if not genomes:
        return None, None, []
    chosen = next((g for g in genomes if g["name"] == genome), genomes[0])
    if contig and contig != "__all__":
        pieces = [c for c in chosen["contigs"] if c["id"] == contig] or list(chosen["contigs"])
    else:
        pieces = list(chosen["contigs"])
    return chosen, contig or "__all__", pieces


def slice_sequence(bakta_dir, genome="", contig="", start=1, end=1, limit=8000):
    chosen, _sel, pieces = _select_genome_pieces(bakta_dir, genome, contig)
    if not chosen or not pieces:
        return {"start": 1, "end": 1, "seq": ""}
    try:
        start = max(int(start), 1)
        end = max(int(end), start)
    except (TypeError, ValueError):
        start, end = 1, 1
    if end - start + 1 > limit:
        end = start + limit - 1
    chunks = []
    offset = 0
    for piece in pieces:
        length = max(int(piece.get("length") or 0), 1)
        p0 = offset + 1
        p1 = offset + length
        a = max(start, p0)
        b = min(end, p1)
        if a <= b:
            seq = piece.get("seq") or ""
            i0 = a - p0
            i1 = b - p0 + 1
            if seq:
                chunks.append(seq[i0:i1])
            else:
                chunks.append("N" * (i1 - i0))
        offset += length
    return {"start": start, "end": end, "seq": "".join(chunks)}
