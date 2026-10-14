"""Parse Torsten Seemann mlst TSV output for the Post-assembly UI."""
from __future__ import annotations

import csv
import os
from collections import Counter
from pathlib import Path


def find_mlst_tsv(out_dir: str | Path) -> Path | None:
    root = Path(out_dir)
    candidates = [
        root / "mlst_out" / "mlst.tsv",
        root / "mlst_out" / "mlst_summary.tsv",
    ]
    for path in candidates:
        if path.is_file() and path.stat().st_size > 0:
            return path
    # Resume / Nextflow publish quirks: search one level deeper.
    nested = root / "mlst_out"
    if nested.is_dir():
        for path in nested.rglob("mlst.tsv"):
            if path.is_file() and path.stat().st_size > 0:
                return path
    return None


def mlst_run_complete(out_dir: str | Path) -> bool:
    return find_mlst_tsv(out_dir) is not None


def _stem(name: str) -> str:
    base = os.path.basename(name or "")
    for ext in (".fasta", ".fa", ".fna", ".FASTA", ".FA", ".FNA"):
        if base.endswith(ext):
            return base[: -len(ext)]
    return base or name


def _dedupe_genomes(rows: list[dict]) -> list[dict]:
    """Keep first row per genome stem (symlink/double-scan can duplicate)."""
    out: list[dict] = []
    seen: set[str] = set()
    for row in rows:
        key = (row.get("genome") or "").strip().lower()
        if not key or key in seen:
            continue
        seen.add(key)
        out.append(row)
    return out


def load_mlst_rows(out_dir: str | Path) -> list[dict]:
    path = find_mlst_tsv(out_dir)
    if not path:
        return []
    rows: list[dict] = []
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        # Full mlst.tsv: FILE SCHEME ST allele...
        # Summary: genome scheme ST
        file_key = "FILE" if "FILE" in fieldnames else ("genome" if "genome" in fieldnames else None)
        scheme_key = "SCHEME" if "SCHEME" in fieldnames else ("scheme" if "scheme" in fieldnames else None)
        st_key = "ST" if "ST" in fieldnames else None
        if not file_key or not scheme_key or not st_key:
            handle.seek(0)
            plain = csv.reader(handle, delimiter="\t")
            header = next(plain, None)
            if not header or len(header) < 3:
                return []
            allele_names = header[3:]
            for cols in plain:
                if len(cols) < 3:
                    continue
                alleles = {
                    allele_names[i]: cols[3 + i]
                    for i in range(min(len(allele_names), len(cols) - 3))
                }
                rows.append({
                    "genome": _stem(cols[0]),
                    "scheme": cols[1] or "-",
                    "ST": cols[2] or "-",
                    "alleles": alleles,
                })
            return _dedupe_genomes(rows)
        allele_keys = [k for k in fieldnames if k not in (file_key, scheme_key, st_key)]
        for raw in reader:
            alleles = {k: (raw.get(k) or "-") for k in allele_keys}
            rows.append({
                "genome": _stem(raw.get(file_key) or ""),
                "scheme": (raw.get(scheme_key) or "-") or "-",
                "ST": (raw.get(st_key) or "-") or "-",
                "alleles": alleles,
            })
    return _dedupe_genomes(rows)


def mlst_display_rows(out_dir: str | Path) -> list[dict]:
    """Flat rows for the DataTable (genome, scheme, ST, allele columns)."""
    out = []
    for row in load_mlst_rows(out_dir):
        flat = {
            "genome": row["genome"],
            "scheme": row["scheme"],
            "ST": row["ST"],
        }
        for gene, allele in (row.get("alleles") or {}).items():
            flat[gene] = allele
        out.append(flat)
    return out


def mlst_chart_payload(rows: list[dict] | None = None, out_dir: str | Path | None = None) -> dict:
    """ST counts + allele matrix for Plotly."""
    if rows is None and out_dir is not None:
        rows = load_mlst_rows(out_dir)
    rows = rows or []
    st_counts = Counter((r.get("ST") or "-") for r in rows)
    scheme_counts = Counter((r.get("scheme") or "-") for r in rows)
    allele_genes: list[str] = []
    seen = set()
    for r in rows:
        for gene in (r.get("alleles") or {}):
            if gene not in seen:
                seen.add(gene)
                allele_genes.append(gene)
    genomes = [r["genome"] for r in rows]
    z = []
    for r in rows:
        alleles = r.get("alleles") or {}
        row_vals = []
        for gene in allele_genes:
            raw = str(alleles.get(gene) or "-")
            # Numeric alleles for colour; non-numeric → 0.
            try:
                row_vals.append(float(raw.split("~")[0].lstrip("~")))
            except ValueError:
                row_vals.append(0.0)
        z.append(row_vals)
    return {
        "st_labels": list(st_counts.keys()),
        "st_counts": [st_counts[k] for k in st_counts],
        "scheme_labels": list(scheme_counts.keys()),
        "scheme_counts": [scheme_counts[k] for k in scheme_counts],
        "genomes": genomes,
        "allele_genes": allele_genes,
        "allele_matrix": z,
        "n_typed": sum(1 for r in rows if (r.get("ST") or "-") not in ("-", "", "0", "-")),
        "n_genomes": len(rows),
    }
