"""Remove Nextflow work directories and other temp files under an output directory."""

from __future__ import annotations

import os
import shutil


TEMP_DIR_NAMES = (
    ".nextflow-work",
    ".nextflow",
    "work",
    ".illumina_filt_work",
)

TEMP_FILE_SUFFIXES = (
    ".partial",
    ".tmp",
    ".temp",
)


def cleanup_run_temps(out_dir: str) -> str:
    """Delete Nextflow work dirs and obvious temp files under out_dir.

    Keeps final results (FASTAs, reports, databases outputs). Returns a
    short status message for the live log.
    """
    root = os.path.abspath(os.path.expanduser(os.path.expandvars(str(out_dir or "").strip())))
    if not root or not os.path.isdir(root):
        return f"BactFlow: cleanup skipped (output directory not found: {out_dir})."

    removed_dirs: list[str] = []
    removed_files = 0

    for name in TEMP_DIR_NAMES:
        path = os.path.join(root, name)
        if os.path.isdir(path):
            try:
                shutil.rmtree(path)
                removed_dirs.append(name)
            except OSError as exc:
                return f"BactFlow: cleanup failed removing {path}: {exc}"

    for dirpath, _dirnames, filenames in os.walk(root):
        # Do not walk into result trees looking for random .tmp that might be tool data;
        # only strip known temp suffixes at the out_dir top level and one level down.
        rel = os.path.relpath(dirpath, root)
        depth = 0 if rel == "." else rel.count(os.sep) + 1
        if depth > 2:
            continue
        for name in filenames:
            low = name.lower()
            if not any(low.endswith(suf) for suf in TEMP_FILE_SUFFIXES):
                continue
            path = os.path.join(dirpath, name)
            try:
                os.remove(path)
                removed_files += 1
            except OSError:
                continue

    # Stale Nextflow log copies often sit next to the pipeline cwd; prefer out_dir.
    for name in (".nextflow.log", ".nextflow.log.1", ".nextflow.log.2"):
        path = os.path.join(root, name)
        if os.path.isfile(path):
            try:
                os.remove(path)
                removed_files += 1
            except OSError:
                pass

    parts = []
    if removed_dirs:
        parts.append("removed " + ", ".join(removed_dirs))
    if removed_files:
        parts.append(f"removed {removed_files} temp file(s)")
    if not parts:
        return f"BactFlow: cleanup complete (nothing to remove under {root})."
    return f"BactFlow: cleanup complete — {'; '.join(parts)} under {root}."
