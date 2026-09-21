#!/bin/env python

from pycirclize import Circos
import numpy as np
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import glob
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqUtils import gc_fraction
import os
import argparse

GBK_EXTS = {".gbk", ".gbff", ".gb"}
GFF_EXTS = {".gff", ".gff3"}
ANNOT_EXTS = GBK_EXTS | GFF_EXTS


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


def _shift_location(loc, offset):
    if loc is None:
        return None
    parts = getattr(loc, "parts", None)
    if parts and len(parts) > 1:
        shifted = [
            FeatureLocation(int(part.start) + offset, int(part.end) + offset, strand=part.strand)
            for part in parts
        ]
        return CompoundLocation(shifted)
    return FeatureLocation(int(loc.start) + offset, int(loc.end) + offset, strand=loc.strand)


def _parse_gff_attrs(col):
    attrs = {}
    for part in (col or "").split(";"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        attrs[key.strip()] = [value.strip()]
    return attrs


def _nice_tick_interval(length):
    """Round tick step aiming for ~12–16 labels along the genome."""
    length = max(int(length), 1)
    target = max(length / 14.0, 1.0)
    nice = (
        5_000, 10_000, 20_000, 25_000, 50_000,
        100_000, 200_000, 250_000, 500_000,
        1_000_000, 2_000_000, 5_000_000, 10_000_000,
    )
    for step in nice:
        if step >= target:
            return step
    return nice[-1]


def _tick_label(value, length):
    value = int(value)
    if length >= 1_000_000:
        return f"{value / 1e6:.2f} Mb" if value else "0 Mb"
    if length >= 10_000:
        return f"{int(round(value / 1e3))} kb" if value else "0 kb"
    return f"{value} bp"


def _tick_positions(length):
    length = max(int(length), 1)
    step = _nice_tick_interval(length)
    positions = list(range(0, length, step))
    if not positions:
        positions = [0]
    near = max(int(length * 0.025), 1)
    trimmed = [p for p in positions if p == 0 or (length - p) > near]
    if trimmed[-1] != length:
        trimmed.append(length)
    return trimmed


class AnnotationSource:
    """All contigs from a Bakta/Prokka folder, concatenated onto one circle."""

    def __init__(self, path):
        self.path = path
        ext = os.path.splitext(path)[1].lower()
        self.name = os.path.splitext(os.path.basename(path))[0]
        self._features = []
        self._seq = ""
        if ext in GBK_EXTS:
            self._load_gbk(path)
        elif ext in GFF_EXTS:
            self._load_gff(path)
        else:
            raise ValueError(f"Unsupported annotation file: {path}")

    def _load_gbk(self, path):
        offset = 0
        seq_parts = []
        for rec in SeqIO.parse(path, "genbank"):
            seq_parts.append(str(rec.seq))
            for feat in rec.features:
                if feat.type in ("source",):
                    continue
                shifted = _shift_location(feat.location, offset)
                if shifted is None:
                    continue
                self._features.append(SeqFeature(shifted, type=feat.type, qualifiers=feat.qualifiers))
            offset += len(rec.seq)
        self._seq = "".join(seq_parts)
        fasta = _companion_fasta(path)
        if (not self._seq.strip("N")) and fasta:
            self._seq = "".join(str(r.seq) for r in SeqIO.parse(fasta, "fasta"))
        if offset == 0 and self._seq:
            offset = len(self._seq)
        self._length = offset or len(self._seq) or 1

    def _load_gff(self, path):
        fasta = _companion_fasta(path)
        contig_order = []
        contig_len = {}
        with open(path, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("##sequence-region"):
                    parts = line.split()
                    if len(parts) >= 4:
                        seqid, start, end = parts[1], int(parts[2]), int(parts[3])
                        contig_order.append(seqid)
                        contig_len[seqid] = max(end - start + 1, 0)
                elif line.startswith("##FASTA"):
                    break

        if fasta:
            fa_recs = [(r.id.split()[0], str(r.seq)) for r in SeqIO.parse(fasta, "fasta")]
            if fa_recs:
                if not contig_order:
                    contig_order = [name for name, _ in fa_recs]
                for name, seq in fa_recs:
                    contig_len[name] = len(seq)
                seq_map = dict(fa_recs)
                self._seq = "".join(seq_map.get(name, "") for name in contig_order)
                extra = [name for name, _ in fa_recs if name not in contig_order]
                if extra:
                    contig_order.extend(extra)
                    self._seq += "".join(seq_map[name] for name in extra)

        offsets = {}
        offset = 0
        for seqid in contig_order:
            offsets[seqid] = offset
            offset += contig_len.get(seqid, 0)

        with open(path, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("##FASTA"):
                    break
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 7:
                    continue
                seqid, ftype, start, end, strand = cols[0], cols[2], int(cols[3]), int(cols[4]), cols[6]
                if seqid not in offsets:
                    offsets[seqid] = offset
                    contig_order.append(seqid)
                    contig_len[seqid] = max(end, contig_len.get(seqid, 0))
                    offset += contig_len[seqid]
                bio_strand = -1 if strand == "-" else 1
                loc = FeatureLocation(offsets[seqid] + start - 1, offsets[seqid] + end, strand=bio_strand)
                attrs = _parse_gff_attrs(cols[8] if len(cols) > 8 else "")
                self._features.append(SeqFeature(loc, type=ftype, qualifiers=attrs))

        self._length = offset or len(self._seq) or 1
        if not self._seq:
            self._seq = "N" * self._length

    @property
    def genome_length(self):
        return int(self._length)

    def extract_features(self, feature_type="CDS", target_strand=None):
        wanted = {t.strip().lower() for t in str(feature_type).split(",") if t.strip()}
        out = []
        for feat in self._features:
            if feat.type.lower() not in wanted:
                continue
            strand = feat.location.strand
            if target_strand is not None and strand != target_strand:
                continue
            out.append(feat)
        return out

    def _window_stats(self, kind):
        seq = (self._seq or "").upper()
        if not seq:
            return np.array([0]), np.array([0.0])
        window_size = 1000 if len(seq) >= 2000 else max(50, len(seq) // 20 or 1)
        pos, vals = [], []
        for i in range(0, max(1, len(seq) - window_size + 1), window_size):
            window = seq[i:i + window_size]
            if kind == "gc":
                vals.append(100.0 * (window.count("G") + window.count("C")) / max(len(window), 1))
            else:
                g, c = window.count("G"), window.count("C")
                denom = g + c
                vals.append(((g - c) / denom) if denom else 0.0)
            pos.append(i + len(window) // 2)
        if not pos:
            pos, vals = [0], [0.0]
        return np.array(pos), np.array(vals)

    def calc_genome_gc_content(self):
        seq = self._seq or ""
        if not seq:
            return 0.0
        return gc_fraction(seq) * 100

    def calc_gc_content(self):
        return self._window_stats("gc")

    def calc_gc_skew(self):
        return self._window_stats("skew")


def collect_annotation_files(root):
    """One annotation file per sample folder, preferring GenBank then GFF."""
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


def annotation_length(path):
    try:
        return AnnotationSource(path).genome_length
    except Exception:
        return 0


# creating the circulize class
class circluar_genome():
    def __init__(self, gbk_path, space, dpi, sector_length=None):
        self.gbk_path = gbk_path
        self.space = space
        self.dpi = dpi
        self.gbk_ref = AnnotationSource(self.gbk_path)
        self.sector_length = int(sector_length or self.gbk_ref.genome_length)
        self.circos = Circos(sectors={self.gbk_ref.name: self.sector_length}, space=self.space)
        self.sector = self.circos.get_sector(self.gbk_ref.name)
        
    def plot_layout(self, text: str, radious = 20, tick_numbers = 10, sector_size = (98,100), fill_color = "lightgrey"):
        
        if text:
            self.circos.text(text, r = radious)
        tic_positions = _tick_positions(self.sector_length)
        outr_track = self.sector.add_track(sector_size)
        outr_track.axis(fc= fill_color)
        outr_track.xticks(
            x=tic_positions,
            labels=[_tick_label(v, self.sector_length) for v in tic_positions],
            label_size=7,
            label_orientation="vertical",
            show_bottom_line=False,
        )
        
    def cds_add(self, gbk = None, sector_size = (95,97), r_pad_ratio= 0.1, name="CDS", color = "black", 
                strand= 1, label = "CDS", label_color = "black", label_size = 3, feature_types=None):
        gbk_obj = AnnotationSource(self.gbk_path if gbk is None else gbk)
        types = feature_types or ([name] if name else ["CDS"])
        feats = []
        for ftype in types:
            feats.extend(gbk_obj.extract_features(feature_type=ftype, target_strand=strand))
        
        if strand ==1:
            f_cds_track = self.sector.add_track(sector_size, r_pad_ratio=r_pad_ratio)
            if feats:
                f_cds_track.genomic_features(feats, fc = color)
            if label is not None:
                self.circos.text( text = label, r = sum(sector_size)/2, size = label_size, color = label_color,  horizontalalignment =  'right')
        elif strand == -1:
            r_cds_track = self.sector.add_track(sector_size, r_pad_ratio=r_pad_ratio)
            if feats:
                r_cds_track.genomic_features(feats, fc = color)
            if label is not None:
                self.circos.text( text = label, r = sum(sector_size)/2, size = label_size, color = label_color,  horizontalalignment =  'right')
        else:
            raise ValueError("Strand must be either +1 or -1")
        
    def gc_add(self, gbk = None, sector_size = (95,97), r_pad_ratio= 0.1, gc_negative_color = "darkgrey", gc_positive_color = "forestgreen", label = "GC content", label_color = "black", label_size = 3):
        gbk_obj = AnnotationSource(self.gbk_path if gbk is None else gbk)

        gc_content_track = self.sector.add_track(sector_size, r_pad_ratio=r_pad_ratio)
        pos_list, gc_contents = gbk_obj.calc_gc_content()
        gc_contents = gc_contents - gbk_obj.calc_genome_gc_content()
        positive_gc_contents = np.where(gc_contents > 0, gc_contents, 0)
        negative_gc_contents = np.where(gc_contents < 0, gc_contents, 0)
        abs_max_gc_content = np.max(np.abs(gc_contents))
        vmin, vmax = -abs_max_gc_content, abs_max_gc_content
        gc_content_track.fill_between(
            pos_list, positive_gc_contents, 0, vmin=vmin, vmax=vmax, color = gc_positive_color
        )
        gc_content_track.fill_between(
            pos_list, negative_gc_contents, 0, vmin=vmin, vmax=vmax, color = gc_negative_color
        )
        self.circos.text( text = label, r = sum(sector_size)/2, size = label_size, color = label_color, horizontalalignment =  'right')
    
    #Adding GC skew
    def gc_skew_add(self, gbk = None, sector_size = (95,97), r_pad_ratio= 0.1,  positive_color = "olive", negative_color = "darkgrey",
        label = "GC skew", label_color = "black", label_size = 3):

        if gbk is None:
           gbk_obj = AnnotationSource(self.gbk_path)
        else:
            gbk_obj = AnnotationSource(gbk)
        gc_skew_track = self.sector.add_track(sector_size, r_pad_ratio=r_pad_ratio)
        pos_list, gc_skews = gbk_obj.calc_gc_skew()
        positive_gc_skew = np.where(gc_skews > 0, gc_skews, 0)
        negative_gc_skew = np.where(gc_skews < 0, gc_skews, 0)
        abs_max_gc_skew = np.max(np.abs(gc_skews))
        vmin, vmax = -abs_max_gc_skew, abs_max_gc_skew
        gc_skew_track.fill_between(
            pos_list, positive_gc_skew, 0, vmin=vmin, vmax=vmax, color = positive_color
        )
        gc_skew_track.fill_between(
            pos_list, negative_gc_skew, 0, vmin=vmin, vmax=vmax, color = negative_color
        )
        self.circos.text( text = label, r = sum(sector_size)/2, size = label_size, color = label_color,   horizontalalignment =  'right')
        
    #adding the variant calls 
    def vc_add(self, vc_df: str, sector_size: (93,95), r_pad_ratio = 0.1, ref_color = "darkmagenta", vc_lab_size = 1.3, vc_color = "olive", ref_catpion = "REF", alt_caption = "ALT", label_size = 3):
        vc_track = self.sector.add_track(sector_size, r_pad_ratio=r_pad_ratio)
        pos_list, ref = vc_df["POS"], vc_df["REF"]
        vc_track.xticks(x = list(pos_list), labels= list(ref), label_size=vc_lab_size, 
                        tick_length=0.60, label_margin = 0.1, label_orientation="horizontal", 
                        text_kws=dict(color = ref_color))
        
        self.circos.text( text = ref_catpion, r =  sector_size[1]+0.9, size = label_size, color = ref_color, horizontalalignment =  'right')
        self.circos.text( text = alt_caption, r = sector_size[0]+0.5, size = label_size, color = vc_color,  horizontalalignment =  'right')

        vc_track = self.sector.add_track(tuple([v + 0.7 for v in sector_size]), r_pad_ratio=r_pad_ratio)
        pos_list, alt = vc_df["POS"], vc_df["ALT"]
        vc_track.xticks(x = list(pos_list), outer=False, labels= list(alt), 
                        label_size=vc_lab_size, tick_length=0, label_orientation="horizontal", 
                        label_margin=0, text_kws=dict(color = vc_color))
        #self.circos.plotfig(dpi=self.dpi)
    
    def save_plot(self, path, dpi = 500):
        self.circos.savefig(path, dpi=dpi)

CIRCLE_OPENING_DEG = 32


def plotter():
    parser = argparse.ArgumentParser(description="Generate a circular genome plot with GC content and skew.")
    parser.add_argument('-d', '--gbk_dir', type= str, required=True, help= "Directory containing Bakta/Prokka folders with .gbk/.gbff or .gff/.gff3 files.")
    parser.add_argument('-o', '--out_dir', type= str, required=True, help= "Output directory.")
    parser.add_argument("--add_gc", default=False, action="store_true", help="Include GC content in the plot.")
    parser.add_argument("--add_skew", action="store_true", default=False, help="Include GC skew in the plot.")
    parser.add_argument("--dpi", type=int, default=300, help="Size of the figure.")
    parser.add_argument("--interval", type=int, default=3, help="Interval for sector size adjustment.")
    parser.add_argument("--figsize", type = int, default=10, help="Figure size (inc)")
    parser.add_argument("--f_color", type=str, default="deepskyblue", help="Forward strand color")
    parser.add_argument("--r_color", type=str, default="#ff7261", help="Reverse strand color")
    parser.add_argument(
        "--feature_types",
        type=str,
        default="cds",
        help="Comma-separated feature types from the gene-type dropdown, e.g. cds,tRNA,rRNA",
    )

    args = parser.parse_args()

    col_forward = args.f_color
    col_reverse = args.r_color
    feature_types = [t.strip() for t in str(args.feature_types).split(",") if t.strip()] or ["cds"]
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)
    out_file = os.path.join(args.out_dir, "circular_plot.png")
    add_gc = args.add_gc
    add_skew = args.add_skew
    dpi = args.dpi
    interval = args.interval

    
    
    gbks = collect_annotation_files(args.gbk_dir)
    lenghts = {}
    for file in gbks:
        lenghts[file] = annotation_length(file)
    if not gbks or not any(lenghts.values()):
        raise SystemExit(f"No .gbk/.gbff/.gff/.gff3 annotation files found under {args.gbk_dir}")
    start_genome = max(lenghts, key=lenghts.get)
    longest_len = max(lenghts.values())

    pl = circluar_genome(start_genome, dpi=dpi, space=CIRCLE_OPENING_DEG, sector_length=longest_len)
    pl.plot_layout(
        text=None 
        )


    ref_name  = os.path.splitext(os.path.basename(start_genome))[0]
    ref_ring = ((97-interval), 97)
    pl.cds_add(
        strand = 1,  
        sector_size= ref_ring, 
        color= col_forward, 
        label=ref_name, 
        label_size = 7,
        feature_types=feature_types,
        )
    pl.cds_add(
        strand = -1, 
        sector_size= ref_ring,  
        color= col_reverse, 
        label=None,
        feature_types=feature_types,
        )


    num_file = len(gbks) - 1

    ranges = []
    start = 97 - interval 

    for _ in range(num_file):
        end = start
        start -=  interval
        ranges.append((start, end))



    #remove the reference
    rm = gbks.index(start_genome)

    gbks.pop(rm)


    for file, index in zip(gbks, range(0, len(gbks))):
        sector_name = os.path.splitext(os.path.basename(file))[0]
        name = sector_name
        
        pl.cds_add(gbk=  file, strand = 1, sector_size=ranges[index], color= col_forward, label = name, label_size= 7, feature_types=feature_types)
        pl.cds_add(gbk=  file, strand = -1, sector_size=ranges[index], color=col_reverse, label=None, feature_types=feature_types)



    #Adding GC

    if add_gc:
        gc_tuple = []
        if ranges:
            start = ranges[-1][-1] - interval
        else:
            start = ref_ring[0] - interval

        n_gc = max(num_file, 1)
        for _ in range(n_gc):
            end = start 
            start -= interval
            gc_tuple.append((start, end)) 

        gc_files = gbks if gbks else [start_genome]
        for file, index in zip(gc_files, range(0, len(gc_files))):
            sector_name = os.path.splitext(os.path.basename(file))[0]
            name = sector_name

            pl.gc_add(label=f"GC content: {name}", sector_size=gc_tuple[index], label_size= 5)

        

    # gc skew
    if add_skew:
        skew_tuple = []
        if add_gc:
            start = gc_tuple[-1][-1] - interval 
        elif ranges:
            start = ranges[-1][-1] - interval
        else:
            start = ref_ring[0] - interval

        skew_files = gbks if gbks else [start_genome]
        for _ in range(len(skew_files)):
            end = start 
            start -= interval
            skew_tuple.append((start, end)) 

        for file, index in zip(skew_files, range(0, len(skew_files))):
            sector_name = os.path.splitext(os.path.basename(file))[0]
            name = sector_name

            pl.gc_skew_add(label=f"GC skew: {name}", sector_size=skew_tuple[index], label_size= 5)



    pl.circos.plotfig(dpi=dpi, figsize=(args.figsize, args.figsize))
    type_label = "/".join(feature_types)
    fwd_label = f"Forward {type_label}"
    rev_label = f"Reverse {type_label}"

    if add_gc and add_skew:
        handles = [
            Patch(color = "white", label =  f"Pangenome of {len(lenghts)} genomes"),
            Patch(color = col_forward, label = fwd_label),
            Patch(color=col_reverse, label = rev_label),
            Line2D([], [], color="forestgreen", label = "Positive GC content", marker="^", ms = 5, ls = "None"),
            Line2D([], [], color="grey", label="Negative GC Content", marker="v", ms=5, ls="None"),
            Line2D([], [], color="olive", label="Positive GC Skew", marker="^", ms=5, ls="None"),
            Line2D([], [], color="grey", label="Negative GC Skew", marker="v", ms=5, ls="None")

        ]
    elif add_gc and not add_skew:
        handles = [
            Patch(color = "white", label =  f"Pangenome of {len(lenghts)} genomes"),
            Patch(color = col_forward, label = fwd_label),
            Patch(color=col_reverse, label = rev_label),
            Line2D([], [], color="forestgreen", label = "Positive GC content", marker="^", ms = 5, ls = "None"),
            Line2D([], [], color="grey", label="Negative GC Content", marker="v", ms=5, ls="None")

        ]
    elif add_skew and not add_gc:
        handles = [
            Patch(color = "white", label =  f"Pangenome of {len(lenghts)} genomes"),
            Patch(color = col_forward, label = fwd_label),
            Patch(color=col_reverse, label = rev_label),
            Line2D([], [], color="olive", label="Positive GC Skew", marker="^", ms=5, ls="None"),
            Line2D([], [], color="grey", label="Negative GC Skew", marker="v", ms=5, ls="None")

        ]
    else:
        handles = [
            Patch(color = "white", label =  f"Pangenome of {len(lenghts)} genomes"),
            Patch(color = col_forward, label = "Forward CDS"),
            Patch(color=col_reverse, label = "Reverse CDS")

        ]

    
    pl.circos.ax.legend(
    handles=handles, 
        bbox_to_anchor = (0.5, 0.5), 
        loc = "center", 
        fontsize = 10
        ).figure.savefig(out_file, dpi = 500)
    print(f"✅ Circular plot saved at {out_file}")
if __name__ == "__main__":
    plotter()