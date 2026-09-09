#!/usr/bin/env python3
"""Parse BactFlow SNP / variant-call VCF folders into table-ready JSON.

Layout expected:
  <out_dir>/snps/<genome>/*.vcf
  <out_dir>/snp/<genome>/*.vcf
  <out_dir>/vcs/<genome>/*.vcf
  <out_dir>/vcf/<genome>/*.vcf

## lines are metadata comments and are skipped.
The column header is the first line that starts with a single '#' and contains CHROM.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
from collections import Counter
from typing import Any


CORE_COLUMNS = ("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")
DEFAULT_MAX_ROWS = 2000


def _open_text(path: str):
    if path.endswith(".gz"):
        import gzip

        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_vcf(path: str, max_rows: int = DEFAULT_MAX_ROWS) -> dict[str, Any]:
    """Parse one VCF into headers + rows. Returns a dict suitable for JSON/UI."""
    result: dict[str, Any] = {
        "path": os.path.abspath(path),
        "genome": os.path.basename(os.path.dirname(path)) or os.path.basename(path),
        "columns": list(CORE_COLUMNS) + ["INFO"],
        "rows": [],
        "n_variants": 0,
        "n_shown": 0,
        "truncated": False,
        "chrom_counts": {},
        "error": None,
    }
    if not os.path.isfile(path):
        result["error"] = f"File not found: {path}"
        return result

    header: list[str] | None = None
    chrom_counts: Counter[str] = Counter()
    rows: list[dict[str, str]] = []
    n_variants = 0

    try:
        with _open_text(path) as handle:
            for raw in handle:
                line = raw.rstrip("\n\r")
                if not line:
                    continue
                if line.startswith("##"):
                    continue
                if line.startswith("#"):
                    # Column header: #CHROM POS ID REF ALT QUAL FILTER INFO ...
                    cols = line.lstrip("#").split("\t")
                    if cols and cols[0].upper() == "CHROM":
                        header = cols
                    continue
                if header is None:
                    continue

                parts = line.split("\t")
                if len(parts) < 8:
                    continue
                n_variants += 1
                chrom = parts[0]
                chrom_counts[chrom] += 1
                if len(rows) >= max_rows:
                    continue
                row = {
                    "CHROM": parts[0],
                    "POS": parts[1],
                    "ID": parts[2],
                    "REF": parts[3],
                    "ALT": parts[4],
                    "QUAL": parts[5],
                    "FILTER": parts[6],
                    "INFO": parts[7][:200] + ("…" if len(parts[7]) > 200 else ""),
                }
                rows.append(row)
    except OSError as exc:
        result["error"] = str(exc)
        return result

    result["n_variants"] = n_variants
    result["n_shown"] = len(rows)
    result["truncated"] = n_variants > len(rows)
    result["rows"] = rows
    result["chrom_counts"] = dict(chrom_counts.most_common(40))
    return result


def find_vcf_files(root: str, basename_only: str | None = None) -> list[str]:
    if not root or not os.path.isdir(root):
        return []
    if basename_only:
        patterns = (
            os.path.join(root, "**", basename_only),
            os.path.join(root, "**", f"{basename_only}.gz"),
        )
    else:
        patterns = (
            os.path.join(root, "**", "*.vcf"),
            os.path.join(root, "**", "*.vcf.gz"),
        )
    found: list[str] = []
    for pattern in patterns:
        found.extend(glob.glob(pattern, recursive=True))
    # Prefer non-gz when both exist; de-dupe
    seen = set()
    out = []
    for path in sorted(found):
        key = path[:-3] if path.endswith(".gz") else path
        if key in seen and path.endswith(".gz"):
            continue
        seen.add(key)
        out.append(path)
    return out


def collect_variant_tables(
    out_dir: str,
    kind: str = "snps",
    max_rows: int = DEFAULT_MAX_ROWS,
) -> dict[str, Any]:
    """Collect parsed VCF tables for SNP or structural/variant-call outputs."""
    kind = (kind or "snps").lower()
    if kind in ("snp", "snps"):
        candidates = ["snps", "snp"]
        label = "SNP"
        # All *.vcf under each genome folder
        name_filter = None
    else:
        candidates = ["vcs", "vcf", "svs", "variants"]
        label = "variant call"
        # Medaka VC: only the primary medaka.vcf (skip annotated/sorted)
        name_filter = "medaka.vcf"

    roots = []
    for name in candidates:
        path = os.path.join(out_dir, name)
        if os.path.isdir(path):
            roots.append(path)

    payload: dict[str, Any] = {
        "exists": False,
        "label": label,
        "roots": roots,
        "genomes": [],
        "error": None,
    }
    if not roots:
        payload["error"] = (
            f"No {label} output folder found under {out_dir} "
            f"(looked for: {', '.join(candidates)})."
        )
        return payload

    genomes = []
    for root in roots:
        for vcf_path in find_vcf_files(root, basename_only=name_filter):
            parsed = parse_vcf(vcf_path, max_rows=max_rows)
            genomes.append(parsed)

    if not genomes:
        wanted = name_filter or "*.vcf"
        payload["error"] = (
            f"{label} folder exists ({', '.join(roots)}) but no {wanted} files were found."
        )
        return payload

    payload["exists"] = True
    payload["genomes"] = genomes
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Summarize BactFlow VCF outputs as JSON tables.")
    parser.add_argument("out_dir", help="BactFlow output directory containing snps/ or vcs/")
    parser.add_argument(
        "--kind",
        choices=("snps", "vcs"),
        default="snps",
        help="Which results folder family to read",
    )
    parser.add_argument("--max-rows", type=int, default=DEFAULT_MAX_ROWS)
    parser.add_argument("-o", "--output", help="Write JSON to this file (default: stdout)")
    args = parser.parse_args(argv)

    payload = collect_variant_tables(args.out_dir, kind=args.kind, max_rows=args.max_rows)
    text = json.dumps(payload, indent=2)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(text + "\n")
    else:
        print(text)
    return 0 if payload.get("exists") else 1


if __name__ == "__main__":
    raise SystemExit(main())
