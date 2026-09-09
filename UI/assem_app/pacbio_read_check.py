#!/usr/bin/env python3
"""Classify PacBio reads as HiFi/CCS vs subreads/CLR.

Does NOT convert FASTQ CLR to HiFi (that requires .subreads.bam + ccs).
For assembly: if reads are subreads/CLR, set Flye mode to pacbio-raw.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
import shutil
import sys
from pathlib import Path

from Bio import SeqIO

SUBREAD_COORD_RE = re.compile(r"/(\d+)_(\d+)(?:\s|$)")
CCS_NAME_RE = re.compile(r"/ccs(?:\s|$)", re.I)
RQ_RE = re.compile(r"(?:^|\s|/)rq:([0-9]*\.?[0-9]+)", re.I)
FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def _open_text(path: Path):
    path = Path(path)
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt", encoding="utf-8", errors="replace")


def _mean_phred(record) -> float | None:
    quals = record.letter_annotations.get("phred_quality")
    if not quals:
        return None
    return float(sum(quals) / len(quals))


def classify_fastq(fastq_path: Path, num_records: int = 10000) -> dict:
    """Classify a PacBio-like FASTQ using headers + mean Phred (+ optional rq)."""
    fastq_path = Path(fastq_path)
    total_qual = 0.0
    n_qual = 0
    record_count = 0
    ccs_header_count = 0
    subread_header_count = 0
    rq_count = 0
    rq_qv_sum = 0.0

    with _open_text(fastq_path) as handle:
        for record in SeqIO.parse(handle, "fastq"):
            record_count += 1
            header = f"{record.id} {record.description}"
            if CCS_NAME_RE.search(header):
                ccs_header_count += 1
            elif SUBREAD_COORD_RE.search(header):
                subread_header_count += 1

            rq_m = RQ_RE.search(header)
            if rq_m:
                rq = min(max(float(rq_m.group(1)), 0.0), 1.0 - 1e-12)
                rq_count += 1
                rq_qv_sum += -10.0 * math.log10(1.0 - rq)

            mean_q = _mean_phred(record)
            if mean_q is not None:
                total_qual += mean_q
                n_qual += 1

            if record_count >= num_records:
                break

    if record_count == 0:
        return {
            "path": str(fastq_path),
            "format": "fastq",
            "read_class": "empty",
            "n_reads_scanned": 0,
            "recommended_flye_mode": None,
            "message": "File appears empty.",
        }

    avg_phred = (total_qual / n_qual) if n_qual else 0.0
    avg_rq_qv = (rq_qv_sum / rq_count) if rq_count else None

    if rq_count >= max(1, record_count // 10) and (avg_rq_qv or 0) >= 20:
        read_class = "hifi"
        reason = "rq tags indicate HiFi/CCS accuracy (QV≥20)"
    elif ccs_header_count > subread_header_count and ccs_header_count > 0:
        read_class = "hifi"
        reason = "headers match PacBio /ccs naming"
    elif subread_header_count > ccs_header_count and subread_header_count > 0:
        read_class = "subread"
        reason = "headers match PacBio subread start_end coordinates"
    elif avg_phred >= 20:
        read_class = "hifi"
        reason = f"mean Phred Q{avg_phred:.1f} ≥ 20 (HiFi-like)"
    elif avg_phred < 15:
        read_class = "subread"
        reason = f"mean Phred Q{avg_phred:.1f} < 15 (CLR/subread-like)"
    else:
        read_class = "unknown"
        reason = f"mixed/ambiguous signal (mean Phred Q{avg_phred:.1f})"

    if read_class == "hifi":
        recommended = "pacbio-hifi"
        user_msg = (
            f"{fastq_path.name}: classified as HiFi/CCS. "
            f"Use Flye mode pacbio-hifi (or pacbio-corr if already corrected). ({reason})"
        )
    elif read_class == "subread":
        recommended = "pacbio-raw"
        user_msg = (
            f"{fastq_path.name}: these look like PacBio subreads/CLR, not HiFi/CCS. "
            f"Set PacBio read type to pacbio-raw for assembly. ({reason})"
        )
    else:
        recommended = "pacbio-raw"
        user_msg = (
            f"{fastq_path.name}: PacBio type is ambiguous. "
            f"Prefer pacbio-raw unless you know these are HiFi/CCS. ({reason})"
        )

    return {
        "path": str(fastq_path),
        "format": "fastq",
        "read_class": read_class,
        "n_reads_scanned": record_count,
        "avg_phred": round(avg_phred, 2),
        "ccs_header_count": ccs_header_count,
        "subread_header_count": subread_header_count,
        "rq_count": rq_count,
        "avg_rq_qv": None if avg_rq_qv is None else round(avg_rq_qv, 2),
        "recommended_flye_mode": recommended,
        "reason": reason,
        "message": user_msg,
    }


def classify_bam(bam_path: Path, num_records: int = 5000) -> dict:
    bam_path = Path(bam_path)
    try:
        import pysam
    except ImportError as exc:
        raise SystemExit("pysam is required to classify BAM files") from exc

    ccs_header_count = 0
    subread_header_count = 0
    rq_count = 0
    n = 0
    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
        for aln in bam.fetch(until_eof=True):
            n += 1
            name = aln.query_name or ""
            if CCS_NAME_RE.search(name) or name.lower().endswith("/ccs"):
                ccs_header_count += 1
            elif SUBREAD_COORD_RE.search(name):
                subread_header_count += 1
            try:
                rq = float(aln.get_tag("rq"))
                rq_count += 1
                if rq >= 0.99:
                    ccs_header_count += 1
            except KeyError:
                pass
            if n >= num_records:
                break

    name_l = bam_path.name.lower()
    if name_l.endswith(".subreads.bam") or subread_header_count > ccs_header_count:
        read_class = "subread"
        reason = "BAM looks like PacBio subreads/CLR"
        recommended = "pacbio-raw"
        user_msg = (
            f"{bam_path.name}: subreads BAM (not HiFi/CCS FASTQ). "
            "Export/convert to HiFi first, or assemble CLR with pacbio-raw. "
            f"({reason})"
        )
    elif ccs_header_count > 0 or rq_count > 0:
        read_class = "hifi"
        reason = "BAM looks like CCS/HiFi reads"
        recommended = "pacbio-hifi"
        user_msg = (
            f"{bam_path.name}: classified as HiFi/CCS BAM. "
            f"Export to FASTQ and use pacbio-hifi. ({reason})"
        )
    else:
        read_class = "unknown"
        reason = "BAM without clear HiFi/subread markers"
        recommended = "pacbio-raw"
        user_msg = (
            f"{bam_path.name}: PacBio BAM type unclear. Prefer pacbio-raw. ({reason})"
        )

    return {
        "path": str(bam_path),
        "format": "bam",
        "read_class": read_class,
        "n_reads_scanned": n,
        "ccs_header_count": ccs_header_count,
        "subread_header_count": subread_header_count,
        "rq_count": rq_count,
        "recommended_flye_mode": recommended,
        "reason": reason,
        "message": user_msg,
    }


def classify_path(path: Path, num_records: int = 10000) -> dict:
    path = Path(path)
    if path.name.lower().endswith(".bam"):
        return classify_bam(path, num_records=min(num_records, 5000))
    return classify_fastq(path, num_records=num_records)


def iter_inputs(path: Path):
    path = Path(path)
    if path.is_file():
        yield path
        return
    if not path.is_dir():
        raise SystemExit(f"Not a file or directory: {path}")
    files = []
    for p in path.iterdir():
        name = p.name.lower()
        if not p.is_file():
            continue
        if name.endswith(".bam") or name.endswith(FASTQ_SUFFIXES):
            files.append(p)
    for p in sorted(files, key=lambda x: x.name.lower()):
        yield p


def consensus_class(results: list[dict]) -> tuple[str, str]:
    classes = {r.get("read_class") for r in results if r.get("read_class") not in (None, "empty")}
    if classes == {"hifi"}:
        return "hifi", "pacbio-hifi"
    if "subread" in classes:
        return "subread", "pacbio-raw"
    if "unknown" in classes:
        return "unknown", "pacbio-raw"
    return "unknown", "pacbio-raw"


def user_facing_summary(results: list[dict]) -> str:
    cls, mode = consensus_class(results)
    lines = [r.get("message") or r.get("reason") or "" for r in results]
    lines = [ln for ln in lines if ln]
    if cls == "subread":
        header = (
            "PacBio check: reads look like SUBREADS/CLR (not HiFi/CCS). "
            f"Set PacBio Flye mode to {mode} before assembly."
        )
    elif cls == "hifi":
        header = (
            "PacBio check: reads look like HiFi/CCS. "
            f"Use Flye mode {mode}."
        )
    else:
        header = (
            "PacBio check: read type is ambiguous. "
            f"Prefer Flye mode {mode} unless you know these are HiFi."
        )
    body = " | ".join(lines[:5])
    return header if not body else f"{header} Details: {body}"


def prepare_passthrough(
    input_path: Path,
    out_dir: Path,
    requested_mode: str,
    num_records: int = 10000,
    fail_on_mode_mismatch: bool = True,
) -> int:
    """Classify inputs, copy FASTQs to out_dir, write messages; enforce mode."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    results = []
    n_fastq = 0
    for path in iter_inputs(input_path):
        info = classify_path(path, num_records=num_records)
        results.append(info)
        print(info.get("message", ""), flush=True)
        name = path.name.lower()
        if name.endswith(".bam"):
            print(
                f"WARNING: skipping BAM {path.name} for Flye "
                "(provide HiFi/CLR FASTQ.gz for assembly).",
                file=sys.stderr,
            )
            continue
        dest = out_dir / path.name
        if path.resolve() != dest.resolve():
            shutil.copy2(path, dest)
        n_fastq += 1

    if n_fastq == 0:
        raise SystemExit(
            "No PacBio FASTQ/FASTQ.gz files found to assemble. "
            "BAM-only inputs are not passed to Flye."
        )

    cls, recommended = consensus_class(results)
    summary_msg = user_facing_summary(results)
    (out_dir / "pacbio_message.txt").write_text(summary_msg + "\n")
    (out_dir / "pacbio_classification.json").write_text(json.dumps(results, indent=2))
    (out_dir / "recommended_flye_mode.txt").write_text(recommended + "\n")
    (out_dir / "read_class.txt").write_text(cls + "\n")
    print(summary_msg, flush=True)

    requested = (requested_mode or "").strip() or "pacbio-hifi"
    if fail_on_mode_mismatch and cls == "subread" and requested == "pacbio-hifi":
        raise SystemExit(
            "ERROR: These PacBio reads are classified as SUBREADS/CLR, not HiFi/CCS. "
            "In the Assembly UI set PacBio read type to pacbio-raw and re-run. "
            f"({summary_msg})"
        )
    if fail_on_mode_mismatch and cls == "hifi" and requested == "pacbio-raw":
        print(
            "WARNING: reads look HiFi/CCS but Flye mode is pacbio-raw; "
            "consider switching to pacbio-hifi.",
            file=sys.stderr,
        )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_cls = sub.add_parser("classify", help="Classify one FASTQ/BAM or a directory")
    p_cls.add_argument("-i", "--input", required=True)
    p_cls.add_argument("-n", "--num-records", type=int, default=10000)
    p_cls.add_argument("--json", action="store_true")

    p_prep = sub.add_parser(
        "prepare",
        help="Classify, copy FASTQs for assembly, enforce recommended Flye mode",
    )
    p_prep.add_argument("-i", "--input", required=True)
    p_prep.add_argument("-o", "--outdir", required=True)
    p_prep.add_argument(
        "--requested-mode",
        default="pacbio-hifi",
        help="UI/Flye mode the user selected (pacbio-raw|pacbio-corr|pacbio-hifi)",
    )
    p_prep.add_argument("-n", "--num-records", type=int, default=10000)
    p_prep.add_argument(
        "--allow-mode-mismatch",
        action="store_true",
        help="Do not exit when subreads are paired with pacbio-hifi",
    )

    args = parser.parse_args(argv)

    if args.cmd == "classify":
        path = Path(args.input)
        if path.is_dir():
            results = [classify_path(p, num_records=args.num_records) for p in iter_inputs(path)]
            payload = {
                "results": results,
                "summary": user_facing_summary(results),
                "consensus_class": consensus_class(results)[0],
                "recommended_flye_mode": consensus_class(results)[1],
            }
            if args.json:
                print(json.dumps(payload, indent=2))
            else:
                print(payload["summary"])
                for r in results:
                    print(f"- {Path(r['path']).name}: {r.get('read_class')} ({r.get('reason')})")
            return 0

        info = classify_path(path, num_records=args.num_records)
        if args.json:
            print(json.dumps(info, indent=2))
        else:
            print(info.get("message") or info.get("reason"))
            for key in (
                "path",
                "read_class",
                "avg_phred",
                "recommended_flye_mode",
                "reason",
            ):
                if key in info and info[key] is not None:
                    print(f"{key}: {info[key]}")
        return 0

    if args.cmd == "prepare":
        return prepare_passthrough(
            Path(args.input),
            Path(args.outdir),
            requested_mode=args.requested_mode,
            num_records=args.num_records,
            fail_on_mode_mismatch=not args.allow_mode_mismatch,
        )

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
