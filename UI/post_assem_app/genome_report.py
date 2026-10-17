"""Deterministic post-assembly report (HTML + PDF) from existing BactFlow outputs."""
from __future__ import annotations

import csv
import html
import os
import re
import shutil
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


WIDE_COL_THRESHOLD = 8
ROWS_PER_PAGE = 28


def _find_logo() -> Path | None:
    """Resolve bactflow logo from UI static paths (host or Docker image layout)."""
    here = Path(__file__).resolve().parent
    candidates = [
        here / "static" / "images" / "logo.png",
        here.parent / "bactflow_app" / "static" / "images" / "logo.png",
        here.parent / "assem_app" / "static" / "images" / "logo.png",
        Path("/app/UI/bactflow_app/static/images/logo.png"),
        Path("/app/UI/post_assem_app/static/images/logo.png"),
    ]
    for path in candidates:
        if path.is_file() and path.stat().st_size > 0:
            return path
    return None


def _stage_logo(assets_dir: Path) -> dict:
    """Copy logo into report assets/; return rel/abs paths or empty dict."""
    src = _find_logo()
    if not src:
        return {}
    assets_dir.mkdir(parents=True, exist_ok=True)
    dest = assets_dir / "bactflow_logo.png"
    shutil.copy2(src, dest)
    return {
        "logo_abs": str(dest),
        "logo_rel": "assets/bactflow_logo.png",
    }


def _read_tsv(path: Path) -> list[dict]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _quast_rows(out_dir: Path) -> list[dict]:
    for name in ("transposed_report.tsv", "report.tsv"):
        path = out_dir / "quast_stat" / name
        rows = _read_tsv(path)
        if rows:
            return rows
    return []


def _checkm_rows(out_dir: Path) -> list[dict]:
    try:
        from checkm_report import find_checkm_files, parse_checkm_lineage_table, read_text_if_nonempty
        lineage_path, _, _ = find_checkm_files(str(out_dir))
        text = read_text_if_nonempty(lineage_path)
        return parse_checkm_lineage_table(text) if text else []
    except Exception:
        return []


def _gtdb_rows(out_dir: Path) -> list[dict]:
    try:
        from gtdb_report import load_gtdb_taxonomy, taxonomy_display_rows
        df = load_gtdb_taxonomy(str(out_dir))
        return taxonomy_display_rows(df) if df is not None and not df.empty else []
    except Exception:
        return []


def _mlst_rows(out_dir: Path) -> list[dict]:
    try:
        from mlst_report import mlst_display_rows
        return mlst_display_rows(str(out_dir))
    except Exception:
        return []


def _plasmid_rows(out_dir: Path) -> tuple[list[dict], list[dict]]:
    try:
        from plasmid_report import plasmid_display_rows, virus_display_rows
        return plasmid_display_rows(str(out_dir)), virus_display_rows(str(out_dir))
    except Exception:
        return [], []


def _bakta_amr_hits(out_dir: Path) -> list[dict]:
    """Collect antibiotic AMRFinderPlus hits from Bakta JSON/TSV if present."""
    bakta = out_dir / "bakta_out"
    if not bakta.is_dir():
        return []
    hits: list[dict] = []
    for tsv in bakta.rglob("*.tsv"):
        name = tsv.name.lower()
        if name.endswith(".hypotheticals.tsv") or name == "gene_count.tsv":
            continue
        try:
            rows = _read_tsv(tsv)
        except Exception:
            continue
        if not rows:
            continue
        keys = {k.lower(): k for k in rows[0].keys()}
        gene_k = keys.get("gene") or keys.get("product")
        type_k = keys.get("type") or keys.get("feature type") or keys.get("feature_type")
        for row in rows:
            typ = str(row.get(type_k) or "").lower() if type_k else ""
            product = str(row.get(keys.get("product") or "") or "").lower()
            gene = str(row.get(gene_k) or "") if gene_k else ""
            blob = " ".join(str(v) for v in row.values()).lower()
            if "amr" in typ or "resistance" in product or "amrfinder" in blob:
                hits.append({
                    "genome": tsv.parent.name or tsv.stem,
                    "gene": gene or row.get(keys.get("product") or "", "-"),
                    "product": row.get(keys.get("product") or "", "-"),
                })
        if len(hits) > 200:
            break
    return hits[:200]


def _collect_figures(out_dir: Path, assets_dir: Path) -> list[dict]:
    """Copy plot assets into report assets/ for HTML; PDF embeds rasters only."""
    assets_dir.mkdir(parents=True, exist_ok=True)
    figures: list[dict] = []
    seen: set[str] = set()

    def add(src: Path, title: str, group: str) -> None:
        if not src.is_file() or src.stat().st_size == 0:
            return
        key = str(src.resolve())
        if key in seen:
            return
        seen.add(key)
        safe = re.sub(r"[^A-Za-z0-9._-]+", "_", f"{group}_{src.name}")[:120]
        dest = assets_dir / safe
        shutil.copy2(src, dest)
        figures.append({
            "title": title,
            "group": group,
            "src": str(src),
            "rel": f"assets/{safe}",
            "abs": str(dest),
            "ext": src.suffix.lower(),
        })

    add(out_dir / "circular_plot.png", "Circular genome map", "circular")

    bakta = out_dir / "bakta_out"
    if bakta.is_dir():
        for png in sorted(bakta.rglob("*.png")):
            add(png, f"Bakta annotation map — {png.stem}", "bakta")

    quast_img = out_dir / "quast_stat" / "basic_stats"
    if quast_img.is_dir():
        for img in sorted(quast_img.iterdir()):
            if img.suffix.lower() in {".png", ".svg", ".jpg", ".jpeg"}:
                add(img, f"QUAST — {img.stem}", "quast")

    for pattern in (
        "*qc*.png",
        "*heatmap*.png",
        "*tree*.png",
        "*checkm*.png",
        "*gtdb*.png",
        "*mlst*.png",
        "requested_genes_*.png",
    ):
        for png in sorted(out_dir.glob(pattern)):
            add(png, png.stem.replace("_", " "), "figures")
        plots = out_dir / "plots"
        if plots.is_dir():
            for png in sorted(plots.glob(pattern)):
                add(png, png.stem.replace("_", " "), "figures")

    for html_path, title in (
        (out_dir / "seqkit_stats.html", "Pre-assembly / seqkit read stats"),
        (out_dir / "quast_stat" / "report.html", "QUAST HTML report"),
        (out_dir / "quast_stat" / "icarus.html", "QUAST Icarus contig viewer"),
    ):
        if html_path.is_file() and html_path.stat().st_size > 0:
            safe = re.sub(r"[^A-Za-z0-9._-]+", "_", html_path.name)[:120]
            dest = assets_dir / safe
            shutil.copy2(html_path, dest)
            figures.append({
                "title": title,
                "group": "html",
                "src": str(html_path),
                "rel": f"assets/{safe}",
                "abs": str(dest),
                "ext": ".html",
            })

    return figures


def _figures_html(figures: list[dict]) -> str:
    imgs = [f for f in figures if f.get("ext") in {".png", ".jpg", ".jpeg", ".svg", ".gif", ".webp"}]
    links = [f for f in figures if f.get("ext") == ".html"]
    if not imgs and not links:
        return ""
    parts = ["<section><h2>Figures</h2>"]
    for fig in imgs:
        parts.append(
            f"<figure><figcaption>{html.escape(fig['title'])}</figcaption>"
            f"<img src='{html.escape(fig['rel'])}' alt='{html.escape(fig['title'])}'></figure>"
        )
    if links:
        parts.append("<h3>Linked HTML reports</h3><ul>")
        for fig in links:
            parts.append(
                f"<li><a href='{html.escape(fig['rel'])}' target='_blank' rel='noopener'>"
                f"{html.escape(fig['title'])}</a></li>"
            )
        parts.append("</ul>")
    parts.append("</section>")
    return "".join(parts)


def _fnum(val, default=None):
    try:
        return float(str(val).replace(",", "").strip())
    except (TypeError, ValueError):
        return default


def _narrative(bundle: dict) -> list[str]:
    paras = []
    n = bundle.get("n_genomes") or 0
    paras.append(
        f"This BactFlow post-assembly report summarises {n} bacterial isolate genome(s) "
        f"analysed on {bundle.get('generated_at', '')}. Values are taken directly from "
        "QUAST, CheckM, GTDB-Tk, MLST (PubMLST), Bakta/AMRFinderPlus, and geNomad where those "
        "steps were run. No external language model was used to invent results."
    )

    quast = bundle.get("quast") or []
    if quast:
        contig_key = next((k for k in quast[0] if re.search(r"#?\s*contigs", k, re.I)), None)
        n50_key = next((k for k in quast[0] if re.search(r"N50", k, re.I)), None)
        len_key = next((k for k in quast[0] if re.search(r"Total length|Genome length", k, re.I)), None)
        contig_vals = [_fnum(r.get(contig_key)) for r in quast if contig_key]
        contig_vals = [v for v in contig_vals if v is not None]
        n50_vals = [_fnum(r.get(n50_key)) for r in quast if n50_key]
        n50_vals = [v for v in n50_vals if v is not None]
        bits = ["Assembly metrics (QUAST)."]
        if contig_vals:
            bits.append(
                f"Contig counts range from {int(min(contig_vals))} to {int(max(contig_vals))} "
                f"(median {int(sorted(contig_vals)[len(contig_vals)//2])})."
            )
            if max(contig_vals) <= 5:
                bits.append("Most assemblies are near-complete (few contigs), consistent with closed or nearly closed chromosomes.")
            elif min(contig_vals) > 50:
                bits.append("Assemblies remain highly fragmented; consider longer reads, higher coverage, or hybrid assembly.")
        if n50_vals:
            bits.append(f"N50 spans {int(min(n50_vals)):,}–{int(max(n50_vals)):,} bp.")
        if len_key:
            lengths = [_fnum(r.get(len_key)) for r in quast]
            lengths = [v for v in lengths if v]
            if lengths:
                bits.append(f"Total assembly length spans {int(min(lengths)):,}–{int(max(lengths)):,} bp.")
        paras.append(" ".join(bits))

    checkm = bundle.get("checkm") or []
    if checkm:
        comp_key = next((k for k in checkm[0] if "Completeness" in k), None)
        cont_key = next((k for k in checkm[0] if "Contamination" in k), None)
        comps = [_fnum(r.get(comp_key)) for r in checkm if comp_key]
        conts = [_fnum(r.get(cont_key)) for r in checkm if cont_key]
        comps = [v for v in comps if v is not None]
        conts = [v for v in conts if v is not None]
        bits = ["Genome quality (CheckM)."]
        if comps:
            good = sum(1 for v in comps if v >= 90)
            bits.append(f"{good}/{len(comps)} genome(s) reach ≥90% completeness.")
        if conts:
            clean = sum(1 for v in conts if v < 5)
            bits.append(f"{clean}/{len(conts)} genome(s) stay below 5% contamination (typical isolate threshold).")
            if any(v >= 10 for v in conts):
                bits.append("One or more genomes show elevated contamination, which can indicate mixed cultures or assembly artifacts.")
        paras.append(" ".join(bits))

    gtdb = bundle.get("gtdb") or []
    if gtdb:
        species_key = next((k for k in gtdb[0] if re.search(r"species|classification|Species", k, re.I)), None)
        species = []
        for r in gtdb:
            if species_key:
                species.append(str(r.get(species_key) or "unclassified"))
            else:
                species.append(str(next(iter(r.values()), "unclassified")))
        top = Counter(species).most_common(5)
        bits = ["Taxonomy (GTDB-Tk)."]
        bits.append(
            "Dominant assignments: "
            + "; ".join(f"{name} (n={count})" for name, count in top)
            + "."
        )
        paras.append(" ".join(bits))

    mlst = bundle.get("mlst") or []
    if mlst:
        typed = [r for r in mlst if str(r.get("ST") or "-") not in ("-", "", "0")]
        schemes = Counter(str(r.get("scheme") or "-") for r in mlst)
        sts = Counter(str(r.get("ST") or "-") for r in typed)
        bits = ["Sequence typing (MLST / PubMLST)."]
        bits.append(f"{len(typed)}/{len(mlst)} genome(s) received a sequence type.")
        if schemes:
            bits.append("Schemes: " + ", ".join(f"{k}×{v}" for k, v in schemes.most_common(6)) + ".")
        if sts:
            bits.append("Most frequent STs: " + ", ".join(f"ST{k} (n={v})" for k, v in sts.most_common(5)) + ".")
        else:
            bits.append("No PubMLST scheme matched these assemblies (common for novel or poorly represented species).")
        paras.append(" ".join(bits))

    amr = bundle.get("amr") or []
    if amr:
        genes = Counter(str(r.get("gene") or r.get("product") or "?") for r in amr)
        bits = [
            f"Antimicrobial resistance (Bakta / AMRFinderPlus): {len(amr)} annotated resistance-associated feature(s).",
            "Top gene labels: " + ", ".join(f"{g} (n={c})" for g, c in genes.most_common(8)) + ".",
            "Review hits in context; metal/stress annotations are not treated as clinical AMR here.",
        ]
        paras.append(" ".join(bits))
    elif bundle.get("bakta_ran"):
        paras.append(
            "Antimicrobial resistance (Bakta / AMRFinderPlus): no antibiotic-resistance features were "
            "extracted from the annotation tables for this report."
        )

    plasmids, viruses = bundle.get("plasmids") or [], bundle.get("viruses") or []
    if plasmids or viruses or bundle.get("genomad_ran"):
        paras.append(
            f"Mobile genetic elements (geNomad): {len(plasmids)} plasmid contig(s) and "
            f"{len(viruses)} virus/provirus contig(s) reported. Empty tables mean the scan ran "
            "but found no high-confidence MGEs."
        )

    paras.append(
        "Interpretation note: completeness, contamination, ST, and AMR calls are only as good as "
        "the assembly and the databases used. Re-run with updated GTDB, Bakta, CheckM, or PubMLST "
        "data when publishing clinical or outbreak conclusions."
    )
    return paras


def collect_report(out_dir: str | Path, assets_dir: Path | None = None) -> dict:
    root = Path(out_dir).resolve()
    quast = _quast_rows(root)
    checkm = _checkm_rows(root)
    gtdb = _gtdb_rows(root)
    mlst = _mlst_rows(root)
    plasmids, viruses = _plasmid_rows(root)
    amr = _bakta_amr_hits(root)
    bakta_ran = (root / "bakta_out").is_dir()
    genomad_ran = (root / "plasmid_out").is_dir()
    if assets_dir is None:
        assets_dir = root / "bactflow_report" / "assets"
    figures = _collect_figures(root, assets_dir)
    logo = _stage_logo(assets_dir)
    genomes = set()
    for rows, key in (
        (quast, None),
        (checkm, None),
        (gtdb, None),
        (mlst, "genome"),
    ):
        for r in rows:
            if key and key in r:
                genomes.add(str(r[key]))
            else:
                if r:
                    genomes.add(str(next(iter(r.values()))))
    bundle = {
        "out_dir": str(root),
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "n_genomes": max(len(genomes), len(quast), len(checkm), len(gtdb), len(mlst)),
        "quast": quast,
        "checkm": checkm,
        "gtdb": gtdb,
        "mlst": mlst,
        "plasmids": plasmids,
        "viruses": viruses,
        "amr": amr,
        "bakta_ran": bakta_ran,
        "genomad_ran": genomad_ran,
        "figures": figures,
        **logo,
    }
    bundle["narrative"] = _narrative(bundle)
    return bundle


def _table_is_wide(rows: list[dict]) -> bool:
    if not rows:
        return False
    cols = list(rows[0].keys())
    if len(cols) >= WIDE_COL_THRESHOLD:
        return True
    # Long cell content also forces landscape
    sample = rows[:5]
    avg = 0
    n = 0
    for row in sample:
        for c in cols:
            avg += len(str(row.get(c, "")))
            n += 1
    return (avg / max(n, 1)) * len(cols) > 180


def _table_html(title: str, rows: list[dict], limit: int = 200) -> str:
    if not rows:
        return f"<section><h2>{html.escape(title)}</h2><p class='muted'>No data.</p></section>"
    cols = list(rows[0].keys())
    wide = _table_is_wide(rows)
    wrap_cls = "table-wrap wide" if wide else "table-wrap"
    head = "".join(f"<th>{html.escape(str(c))}</th>" for c in cols)
    body = []
    for row in rows[:limit]:
        body.append("<tr>" + "".join(f"<td>{html.escape(str(row.get(c, '')))}</td>" for c in cols) + "</tr>")
    note = ""
    if len(rows) > limit:
        note = f"<p class='muted'>Showing first {limit} of {len(rows)} rows.</p>"
    return (
        f"<section>"
        f"<h2>{html.escape(title)}</h2>{note}"
        f"<div class='{wrap_cls}'><table><thead><tr>{head}</tr></thead>"
        f"<tbody>{''.join(body)}</tbody></table></div></section>"
    )




def write_html_report(bundle: dict, path: Path) -> Path:
    paras = "".join(f"<p>{html.escape(p)}</p>" for p in bundle.get("narrative") or [])
    logo_rel = bundle.get("logo_rel") or ""
    logo_html = ""
    if logo_rel:
        logo_html = (
            f"<img class='cover-logo' src='{html.escape(logo_rel)}' "
            f"alt='' width='48' height='48'>"
        )
    cover = f"""
<section class="cover">
  <div class="cover-brand">
    {logo_html}
    <p class="cover-kicker">BactFlow</p>
  </div>
  <h1>Genome report</h1>
  <p class="meta">{html.escape(bundle.get('generated_at', ''))}<br>
  {html.escape(bundle.get('out_dir', ''))}</p>
</section>
"""
    fig_html = _figures_html(bundle.get("figures") or [])
    sections = [
        cover,
        f"<section><h2>Summary</h2>{paras}</section>",
        fig_html,
        _table_html("QUAST assembly metrics", bundle.get("quast") or []),
        _table_html("CheckM completeness", bundle.get("checkm") or []),
        _table_html("GTDB-Tk taxonomy", bundle.get("gtdb") or []),
        _table_html("MLST (PubMLST)", bundle.get("mlst") or []),
        _table_html("AMR features (Bakta / AMRFinderPlus)", bundle.get("amr") or []),
        _table_html("Plasmids (geNomad)", bundle.get("plasmids") or []),
        _table_html("Viruses / proviruses (geNomad)", bundle.get("viruses") or []),
    ]
    doc = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<title>BactFlow genome report</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Source+Serif+4:opsz,wght@8..60,400;8..60,600&display=swap" rel="stylesheet">
<style>
body{{font-family:"Source Serif 4",Georgia,"Iowan Old Style",serif;max-width:1100px;margin:2rem auto;padding:0;color:#1f2937;line-height:1.45;background:#fafafa}}
h1{{color:#0f3d2e;font-size:1.9rem;margin:.2rem 0 .8rem}} h2{{color:#14532d;border-bottom:1px solid #d1d5db;padding-bottom:.3rem;margin-top:2rem}}
h3{{color:#166534;font-size:1.05rem}}
.meta,.muted{{color:#6b7280;font-size:.95rem}}
.cover{{padding:2rem 10px 1.5rem;border-bottom:1px solid #d1d5db;margin:0 0 1.5rem;
background:radial-gradient(120% 80% at 10% 0%,#ecfdf5 0%,transparent 55%),radial-gradient(90% 60% at 100% 20%,#f0fdf4 0%,transparent 50%),#fafafa}}
.cover-brand{{display:flex;align-items:center;justify-content:flex-start;gap:5px;margin:0 0 .6rem;padding:0}}
.cover-logo{{flex:0 0 auto;width:48px;height:48px;object-fit:contain;border:0;background:transparent;display:block}}
.cover-kicker{{letter-spacing:.28em;text-transform:uppercase;font-size:.78rem;color:#166534;margin:0;font-weight:600;line-height:1}}
.cover h1{{font-size:clamp(2.2rem,5vw,3.2rem);font-weight:600;letter-spacing:-.02em;margin:0.2rem 0 0.8rem;text-align:center}}
section:not(.cover),footer{{padding-left:1rem;padding-right:1rem}}
.table-wrap{{overflow-x:auto;margin:1rem 0;-webkit-overflow-scrolling:touch}}
.table-wrap.wide table{{min-width:960px}}
table{{border-collapse:collapse;width:100%;font-size:.82rem}}
th,td{{border:1px solid #d1d5db;padding:.35rem .5rem;text-align:left;vertical-align:top}}
th{{background:#ecfdf5}}
figure{{margin:1.25rem 0;padding:0;page-break-inside:avoid}}
figcaption{{font-size:.95rem;color:#374151;margin-bottom:.4rem;font-weight:600}}
img{{max-width:100%;height:auto;border:1px solid #e5e7eb;background:#fff}}
footer{{margin:2.5rem 0 1rem;padding-top:1rem;border-top:1px solid #e5e7eb}}
@media print{{
  body{{margin:0;background:#fff;max-width:none}}
  .cover{{page-break-after:avoid;background:#fff;padding-left:10px;padding-right:10px}}
  .cover-logo{{border:0;background:transparent}}
  .table-wrap{{overflow:visible}}
  figure{{page-break-inside:avoid}}
}}
</style></head><body>
{''.join(sections)}
<footer><p class='muted'>Generated by BactFlow. Cite the underlying tools (QUAST, CheckM, GTDB-Tk, mlst, Bakta, geNomad) when publishing. Open this HTML next to its <code>assets/</code> folder so figures load.</p></footer>
</body></html>"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(doc, encoding="utf-8")
    return path


def write_pdf_report(bundle: dict, path: Path) -> Path:
    try:
        from reportlab.lib import colors
        from reportlab.lib.pagesizes import A4, landscape
        from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
        from reportlab.lib.units import cm
        from reportlab.platypus import (
            BaseDocTemplate,
            Frame,
            Image,
            KeepTogether,
            NextPageTemplate,
            PageBreak,
            PageTemplate,
            Paragraph,
            Spacer,
            Table,
            TableStyle,
        )
    except ImportError as exc:
        raise RuntimeError(
            "reportlab is required for PDF reports. Install with: pip install reportlab"
        ) from exc

    path.parent.mkdir(parents=True, exist_ok=True)

    # CSS px → PDF points (96 CSS px per inch).
    px = 72.0 / 96.0
    left_brand = 10 * px  # logo + BactFlow sit 10px from the page edge
    left = right = top = bottom = 1.4 * cm
    portrait_w, portrait_h = A4
    land_w, land_h = landscape(A4)
    page_w = land_w - left - right
    page_h = land_h - top - bottom

    doc = BaseDocTemplate(
        str(path),
        pagesize=A4,
        leftMargin=left_brand,
        rightMargin=right,
        topMargin=top,
        bottomMargin=bottom,
        title="BactFlow genome report",
    )
    # Portrait uses a tight left margin so the brand hugs the page edge.
    frame_p = Frame(
        left_brand,
        bottom,
        portrait_w - left_brand - right,
        portrait_h - top - bottom,
        id="portrait",
        leftPadding=0,
        rightPadding=0,
        topPadding=0,
        bottomPadding=0,
    )
    frame_l = Frame(left, bottom, page_w, page_h, id="landscape")
    doc.addPageTemplates([
        PageTemplate(id="portrait", frames=[frame_p], pagesize=A4),
        PageTemplate(id="landscape", frames=[frame_l], pagesize=landscape(A4)),
    ])

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(
        name="BFTitle",
        parent=styles["Title"],
        textColor=colors.HexColor("#0f3d2e"),
        fontSize=26,
        spaceAfter=6,
        alignment=1,  # CENTER
        leftIndent=0,
    ))
    styles.add(ParagraphStyle(
        name="BFKicker",
        parent=styles["Normal"],
        textColor=colors.HexColor("#166534"),
        fontSize=9,
        letterSpacing=2,
        spaceAfter=10,
        alignment=0,
        leftIndent=0,
    ))
    styles.add(ParagraphStyle(name="BFH2", parent=styles["Heading2"], textColor=colors.HexColor("#14532d"), spaceBefore=8, spaceAfter=6, alignment=0))
    styles.add(ParagraphStyle(name="BFH3", parent=styles["Heading3"], textColor=colors.HexColor("#166534"), fontSize=11, spaceBefore=4, spaceAfter=4))
    styles.add(ParagraphStyle(name="BFBody", parent=styles["BodyText"], leading=14, spaceAfter=8))
    styles.add(ParagraphStyle(name="BFMuted", parent=styles["BodyText"], textColor=colors.grey, fontSize=9, alignment=0))
    styles.add(ParagraphStyle(name="BFCell", parent=styles["BodyText"], fontSize=7, leading=9))

    # Portrait: title + summary only. Everything after switches to landscape.
    # Logo immediately left of "BACTFLOW" (5px gap), brand flush left (10px page margin).
    brand_row: list = []
    logo_abs = bundle.get("logo_abs") or ""
    if logo_abs and Path(logo_abs).is_file():
        logo_img = Image(logo_abs)
        max_logo = 1.15 * cm
        iw, ih = float(logo_img.imageWidth or 1), float(logo_img.imageHeight or 1)
        scale = min(max_logo / iw, max_logo / ih, 1.0)
        logo_img.drawWidth = iw * scale
        logo_img.drawHeight = ih * scale
        gap = 5 * px
        brand = Table(
            [[logo_img, Paragraph("BACTFLOW", styles["BFKicker"])]],
            colWidths=[logo_img.drawWidth, portrait_w - left_brand - right - logo_img.drawWidth - gap],
            hAlign="LEFT",
        )
        brand.setStyle(TableStyle([
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("ALIGN", (0, 0), (0, 0), "LEFT"),
            ("ALIGN", (1, 0), (1, 0), "LEFT"),
            ("LEFTPADDING", (0, 0), (0, 0), 0),
            ("RIGHTPADDING", (0, 0), (0, 0), gap),
            ("LEFTPADDING", (1, 0), (1, 0), 0),
            ("RIGHTPADDING", (1, 0), (1, 0), 0),
            ("TOPPADDING", (0, 0), (-1, -1), 0),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
        ]))
        brand_row = [brand]
    else:
        brand_row = [Paragraph("BACTFLOW", styles["BFKicker"])]

    story: list = brand_row + [
        Paragraph("Genome report", styles["BFTitle"]),
        Paragraph(
            f"{html.escape(bundle.get('generated_at', ''))}<br/>{html.escape(bundle.get('out_dir', ''))}",
            styles["BFMuted"],
        ),
        Spacer(1, 0.5 * cm),
        Paragraph("Summary", styles["BFH2"]),
    ]
    for para in bundle.get("narrative") or []:
        story.append(Paragraph(html.escape(para), styles["BFBody"]))

    story.append(NextPageTemplate("landscape"))
    story.append(PageBreak())

    # Raster figures: scale so caption + image fit one landscape page (no orphan caption pages).
    imgs = [
        f for f in (bundle.get("figures") or [])
        if f.get("ext") in {".png", ".jpg", ".jpeg", ".gif", ".webp"}
        and Path(f.get("abs") or "").is_file()
        and Path(f["abs"]).stat().st_size > 0
    ]
    for idx, fig in enumerate(imgs):
        try:
            probe = Image(fig["abs"])
            iw, ih = float(probe.imageWidth), float(probe.imageHeight)
            if iw <= 0 or ih <= 0:
                continue
            # Keep caption + image on one landscape page. Headings need real vertical
            # budget or KeepTogether splits and leaves an empty-looking caption page.
            reserve = 3.6 * cm if idx == 0 else 2.0 * cm
            max_w = page_w * 0.98
            max_h = max(page_h - reserve, 3.5 * cm) * 0.95
            scale = min(max_w / iw, max_h / ih, 1.0)
            draw_w, draw_h = iw * scale, ih * scale
            # Hard cap: never exceed ~78% of frame height.
            if draw_h > page_h * 0.78:
                s2 = (page_h * 0.78) / ih
                draw_w, draw_h = iw * s2, ih * s2
            if draw_w > page_w:
                s3 = page_w / iw
                draw_w, draw_h = iw * s3, ih * s3
            img = Image(fig["abs"], width=draw_w, height=draw_h)
            img.hAlign = "CENTER"
            block = []
            if idx == 0:
                block.append(Paragraph("Figures", styles["BFH2"]))
            block.append(Paragraph(html.escape(fig["title"]), styles["BFH3"]))
            block.append(img)
            story.append(KeepTogether(block))
            if idx < len(imgs) - 1:
                story.append(Spacer(1, 0.25 * cm))
        except Exception:
            continue

    def _styled_table(data: list[list], width: float) -> Table:
        ncols = max(len(data[0]), 1)
        col_w = width / ncols
        wrapped = []
        for row in data:
            wrapped.append([
                Paragraph(html.escape(str(cell)), styles["BFCell"]) for cell in row
            ])
        table = Table(wrapped, colWidths=[col_w] * ncols, repeatRows=1)
        table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#d1fae5")),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE", (0, 0), (-1, -1), 7),
            ("GRID", (0, 0), (-1, -1), 0.25, colors.grey),
            ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#fafafa")]),
            ("LEFTPADDING", (0, 0), (-1, -1), 2),
            ("RIGHTPADDING", (0, 0), (-1, -1), 2),
            ("TOPPADDING", (0, 0), (-1, -1), 2),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
        ]))
        return table

    def add_table(title: str, rows: list[dict]) -> None:
        # Skip empty tables so they do not waste landscape pages.
        if not rows:
            return
        story.append(Paragraph(html.escape(title), styles["BFH2"]))
        cols = list(rows[0].keys())
        header = [str(c) for c in cols]
        for start in range(0, len(rows), ROWS_PER_PAGE):
            chunk = rows[start:start + ROWS_PER_PAGE]
            data = [header] + [[str(row.get(c, "")) for c in cols] for row in chunk]
            story.append(_styled_table(data, page_w))
            more = start + ROWS_PER_PAGE < len(rows)
            if more:
                story.append(Paragraph(
                    f"Rows {start + 1}–{start + len(chunk)} of {len(rows)} (continued…)",
                    styles["BFMuted"],
                ))
                story.append(PageBreak())
                story.append(Paragraph(html.escape(title) + " (continued)", styles["BFH2"]))

    add_table("QUAST assembly metrics", bundle.get("quast") or [])
    add_table("CheckM completeness", bundle.get("checkm") or [])
    add_table("GTDB-Tk taxonomy", bundle.get("gtdb") or [])
    add_table("MLST (PubMLST)", bundle.get("mlst") or [])
    add_table("AMR features (Bakta / AMRFinderPlus)", bundle.get("amr") or [])
    add_table("Plasmids (geNomad)", bundle.get("plasmids") or [])
    add_table("Viruses / proviruses (geNomad)", bundle.get("viruses") or [])

    story.append(Spacer(1, 0.4 * cm))
    story.append(Paragraph(
        "Generated by BactFlow. Cite QUAST, CheckM, GTDB-Tk, mlst, Bakta, and geNomad when publishing. "
        "An HTML report with the same tables and figures is saved alongside this PDF.",
        styles["BFMuted"],
    ))
    doc.build(story)
    return path


def build_reports(out_dir: str | Path) -> dict:
    root = Path(out_dir).resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"Output directory does not exist: {root}")
    report_dir = root / "bactflow_report"
    assets_dir = report_dir / "assets"
    if assets_dir.exists():
        shutil.rmtree(assets_dir)
    assets_dir.mkdir(parents=True, exist_ok=True)
    bundle = collect_report(root, assets_dir=assets_dir)
    html_path = write_html_report(bundle, report_dir / "bactflow_genome_report.html")
    pdf_path = report_dir / "bactflow_genome_report.pdf"
    write_pdf_report(bundle, pdf_path)
    n_img = len([
        f for f in (bundle.get("figures") or [])
        if f.get("ext") in {".png", ".jpg", ".jpeg", ".gif", ".webp"}
    ])
    return {
        "ok": True,
        "out_dir": str(root),
        "report_dir": str(report_dir),
        "html": str(html_path),
        "pdf": str(pdf_path),
        "n_genomes": bundle.get("n_genomes") or 0,
        "n_figures": n_img,
        "sections": {
            "quast": len(bundle.get("quast") or []),
            "checkm": len(bundle.get("checkm") or []),
            "gtdb": len(bundle.get("gtdb") or []),
            "mlst": len(bundle.get("mlst") or []),
            "amr": len(bundle.get("amr") or []),
            "plasmids": len(bundle.get("plasmids") or []),
            "viruses": len(bundle.get("viruses") or []),
            "figures": n_img,
        },
        "narrative": bundle.get("narrative") or [],
    }
