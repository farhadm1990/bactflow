#!/usr/bin/env python3
"""Pool per-sample FASTQ folders (ONT barcodes or PacBio movie/subread dirs).

Layouts handled:
  parent/
    TL110_fastq/*.subreads.fastq   -> pooled/TL110.fastq
    TL19_fastq/*.subreads.fastq    -> pooled/TL19.fastq
  TL110_fastq/*.subreads.fastq     -> pooled/TL110.fastq  (pointed at one sample)
  parent/*.fastq.gz                -> copied into pooled/ as-is (already one file per sample)

Empty files are skipped. Gzipped inputs are decompressed into uncompressed .fastq.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

class PoolReadsError(Exception):
    """User-facing pooling error (missing dir, no FASTQs, etc.)."""


FASTQ_EXTS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
SKIP_DIRS = {"pooled", "illumina_filt"}
PACBIO_MOVIE_RE = re.compile(r"^m\d{6}_", re.I)
SAMPLE_SUFFIXES = (
    "_fastq_pass",
    "_fastq",
    "_fastqs",
    "_reads",
    "_pass",
)


def _lower(name: str) -> str:
    return name.lower()


def is_fastq_name(name: str, extension: str | None = None) -> bool:
    n = _lower(name)
    if n.endswith(".bam"):
        return False
    ext = _normalize_ext(extension)
    if ext and ext not in ("auto",):
        return n.endswith(ext)
    return n.endswith(FASTQ_EXTS)


def _normalize_ext(extension: str | None) -> str:
    ext = (extension or "auto").strip().lower()
    if ext in ("", "auto", "detect", "*"):
        return "auto"
    if not ext.startswith("."):
        ext = "." + ext
    return ext


def sample_name_from_dir(name: str) -> str:
    n = name
    lower = n.lower()
    for suffix in SAMPLE_SUFFIXES:
        if lower.endswith(suffix):
            n = n[: -len(suffix)]
            break
    return n or name


def human_size(num_bytes: int) -> str:
    value = float(num_bytes)
    for unit in ("B", "K", "M", "G", "T"):
        if value < 1024.0 or unit == "T":
            if unit == "B":
                return f"{int(value)}{unit}"
            return f"{value:.1f}{unit}"
        value /= 1024.0
    return f"{num_bytes}B"


def list_read_files(directory, recurse: bool = False) -> list[str]:
    if not directory or not os.path.isdir(directory):
        return []
    files = []
    try:
        names = sorted(os.listdir(directory))
    except OSError:
        return []
    for name in names:
        path = os.path.join(directory, name)
        if os.path.isfile(path) and is_fastq_name(name):
            files.append(path)
        elif recurse and os.path.isdir(path) and name not in SKIP_DIRS:
            files.extend(list_read_files(path, recurse=False))
    return files


def _iter_candidate_files(root: Path) -> list[Path]:
    files = []
    if not root.is_dir():
        return files
    for child in root.iterdir():
        if child.is_file():
            files.append(child)
        elif child.is_dir() and child.name not in SKIP_DIRS:
            for nested in child.iterdir():
                if nested.is_file():
                    files.append(nested)
    return files


def sniff_extension(fastq_dir: str | Path) -> str:
    counts = {ext: 0 for ext in FASTQ_EXTS}
    subreads = 0
    root = Path(fastq_dir)
    for path in _iter_candidate_files(root):
        name = _lower(path.name)
        if "subreads" in name:
            subreads += 1
        for ext in FASTQ_EXTS:
            if name.endswith(ext):
                counts[ext] += 1
                break
    if not any(counts.values()):
        return ".fastq.gz"
    if counts[".fastq"] and counts[".fastq"] >= counts[".fastq.gz"]:
        return ".subreads.fastq" if subreads else ".fastq"
    if counts[".fastq.gz"]:
        return ".subreads.fastq.gz" if subreads and counts[".fastq"] == 0 else ".fastq.gz"
    if counts[".fq.gz"]:
        return ".fq.gz"
    return ".fq"


def _nonempty_fastqs(directory: Path, extension: str) -> list[Path]:
    files = []
    if not directory.is_dir():
        return files
    for path in sorted(directory.iterdir()):
        if not path.is_file():
            continue
        if not is_fastq_name(path.name, extension):
            continue
        try:
            if path.stat().st_size <= 0:
                continue
        except OSError:
            continue
        files.append(path)
    return files


def _looks_like_movie_chunks(files: list[Path]) -> bool:
    if len(files) < 2:
        return False
    subread_n = sum(1 for p in files if "subreads" in _lower(p.name))
    movie_n = sum(1 for p in files if PACBIO_MOVIE_RE.match(p.name))
    return subread_n >= max(2, len(files) // 2) or movie_n >= max(2, len(files) // 2)


def collect_sample_groups(fastq_dir: str | Path, extension: str = "auto") -> list[tuple[str, list[Path]]]:
    root = Path(fastq_dir)
    ext = _normalize_ext(extension)
    nested = []
    for child in sorted(root.iterdir()):
        if not child.is_dir() or child.name in SKIP_DIRS:
            continue
        files = _nonempty_fastqs(child, ext)
        if files:
            nested.append((sample_name_from_dir(child.name), files))
    if nested:
        return nested

    top = _nonempty_fastqs(root, ext)
    if not top:
        return []
    if _looks_like_movie_chunks(top):
        return [(sample_name_from_dir(root.name), top)]
    return [(Path(p.name).name, [p]) for p in top]


def inspect_layout(fastq_dir: str | Path, extension: str = "auto") -> dict:
    root = Path(fastq_dir) if fastq_dir else Path()
    empty = {
        "directory": str(root) if fastq_dir else "",
        "exists": False,
        "extension": ".fastq.gz",
        "nested_samples": [],
        "top_level_files": [],
        "pooled_files": [],
        "looks_like_subreads": False,
        "needs_concat": False,
        "movie_chunk_dir": False,
        "note": "",
    }
    if not fastq_dir or not root.is_dir():
        return empty

    sniffed = sniff_extension(root)
    nested_names = []
    nested_file_count = 0
    for child in sorted(root.iterdir()):
        if not child.is_dir() or child.name in SKIP_DIRS:
            continue
        files = _nonempty_fastqs(child, "auto")
        if files:
            nested_names.append(child.name)
            nested_file_count += len(files)

    top = _nonempty_fastqs(root, "auto")
    pooled = _nonempty_fastqs(root / "pooled", "auto")
    groups = collect_sample_groups(root, "auto")
    movie_dir = bool(top) and _looks_like_movie_chunks(top) and not nested_names
    looks_sub = any("subreads" in _lower(p.name) for p in _iter_candidate_files(root))
    needs_concat = bool(nested_names) or movie_dir
    note = ""
    if nested_names:
        note = (
            f"Found {len(nested_names)} sample folder(s) "
            f"({', '.join(nested_names[:6])}{'…' if len(nested_names) > 6 else ''}) "
            f"with {nested_file_count} FASTQ file(s). Concatenate to pool each sample."
        )
    elif movie_dir:
        note = (
            f"This folder looks like PacBio movie/subread chunks ({len(top)} files). "
            "Concatenate to pool them into one FASTQ before assembly."
        )

    return {
        "directory": str(root),
        "exists": True,
        "extension": sniffed,
        "nested_samples": nested_names,
        "top_level_files": [p.name for p in top],
        "pooled_files": [p.name for p in pooled],
        "looks_like_subreads": looks_sub,
        "needs_concat": needs_concat,
        "movie_chunk_dir": movie_dir,
        "sample_count": len(groups) if needs_concat else len(top) or len(pooled),
        "note": note,
    }


def _concat_paths(sources: list[Path], dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".partial")
    try:
        with open(tmp, "wb") as out:
            for src in sources:
                if _lower(src.name).endswith(".gz"):
                    with gzip.open(src, "rb") as handle:
                        shutil.copyfileobj(handle, out)
                else:
                    with open(src, "rb") as handle:
                        shutil.copyfileobj(handle, out)
        tmp.replace(dest)
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise


def _up_to_date(dest: Path, sources: list[Path]) -> bool:
    if not dest.is_file() or dest.stat().st_size <= 0:
        return False
    try:
        dest_mtime = dest.stat().st_mtime
        return dest_mtime >= max(p.stat().st_mtime for p in sources)
    except OSError:
        return False


def _safe_sample_filename(sample: str) -> str:
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", sample).strip("._") or "sample"
    if name.lower().endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
        for ext in FASTQ_EXTS:
            if name.lower().endswith(ext):
                name = name[: -len(ext)]
                break
        name = name.rstrip(".")
    if name.lower().endswith(".subreads"):
        name = name[: -len(".subreads")]
    return f"{name}.fastq"


def pool_fastq_dir(fastq_dir: str | Path, extension: str = "auto", cpus: int = 1) -> list[Path]:
    root = Path(fastq_dir).resolve()
    if not root.is_dir():
        raise PoolReadsError(f"FASTQ directory does not exist: {fastq_dir}")

    ext = _normalize_ext(extension)
    if ext == "auto":
        ext = sniff_extension(root)

    groups = collect_sample_groups(root, ext)
    if not groups:
        groups = collect_sample_groups(root, "auto")
    if not groups:
        pooled_existing = _nonempty_fastqs(root / "pooled", "auto")
        if pooled_existing:
            print(f"Nothing new to concatenate; using {len(pooled_existing)} file(s) in {root / 'pooled'}")
            return pooled_existing
        raise PoolReadsError(
            f"No FASTQ files found to concatenate in {root}. "
            "Expected sample subfolders (e.g. TL110_fastq/*.subreads.fastq) "
            "or FASTQ files in the directory."
        )

    pooled_dir = root / "pooled"
    pooled_dir.mkdir(parents=True, exist_ok=True)
    jobs = []
    for sample, sources in groups:
        dest = pooled_dir / _safe_sample_filename(sample)
        if _up_to_date(dest, sources):
            print(f"Keeping existing {dest.name} ({len(sources)} source file(s))")
            jobs.append(("skip", dest, sources))
        else:
            jobs.append(("write", dest, sources))

    workers = max(1, min(int(cpus or 1), sum(1 for kind, _, _ in jobs if kind == "write") or 1))

    def _run(job):
        kind, dest, sources = job
        if kind == "write":
            print(f"Pooling {len(sources)} file(s) -> {dest.name}")
            _concat_paths(sources, dest)
        return dest

    written = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(_run, job) for job in jobs]
        for fut in as_completed(futures):
            written.append(fut.result())

    written = sorted(written, key=lambda p: p.name.lower())
    ready = root / "concatenated_fq_are_ready"
    ready.write_text("ok\n")
    print(f"Pooled {len(written)} FASTQ file(s) in {pooled_dir}")
    return written


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-g", "--fastq-dir", required=True, help="FASTQ directory")
    parser.add_argument("-c", "--cpus", type=int, default=1)
    parser.add_argument("-e", "--extension", default="auto")
    parser.add_argument("-o", "--output-dir", default="", help="Ignored; pooled/ is written next to the inputs")
    parser.add_argument("--inspect", action="store_true", help="Print layout JSON and exit")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    if args.inspect:
        import json

        print(json.dumps(inspect_layout(args.fastq_dir, args.extension), indent=2))
        return 0
    try:
        pool_fastq_dir(args.fastq_dir, extension=args.extension, cpus=args.cpus)
    except PoolReadsError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
