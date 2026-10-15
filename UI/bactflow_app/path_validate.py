"""Shared directory preflight for QC, Assembly, and Post-assembly.

Typed paths must exist (or be creatable for outputs). Content-bearing inputs
must also contain the expected files. Typos are suggested, never auto-used.
"""

from __future__ import annotations

import os
from typing import Any


FASTQ_EXTS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
FASTA_EXTS = (".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz")


def _truthy(value, default=False) -> bool:
    if value is None or str(value).strip() == "":
        return default
    return str(value).strip().lower() in ("true", "1", "yes", "on")


def normalize_path(raw: str | None) -> str:
    text = str(raw or "").strip().strip('"').strip("'")
    if not text:
        return ""
    return os.path.abspath(os.path.expanduser(os.path.expandvars(text)))


def _result(ok: bool, message: str, *, path: str = "", **extra: Any) -> dict:
    payload = {"ok": ok, "message": message, "path": path or "", "directory": path or ""}
    payload.update(extra)
    return payload


def suggest_dir(raw: str, *, kind: str = "fastq") -> str:
    """Suggest a likely typo fix without rewriting the typed path."""
    text = str(raw or "").strip().strip('"').strip("'")
    if not text:
        return ""
    path = normalize_path(text)
    if os.path.isdir(path):
        return ""
    candidates = []
    if kind == "fastq":
        try:
            from pool_reads import resolve_fastq_dir
        except ImportError:
            resolve_fastq_dir = None
        if resolve_fastq_dir:
            try:
                alt, note = resolve_fastq_dir(text, allow_alias=True)
                if note and alt:
                    return str(alt)
            except Exception:
                pass
        for old, new in (("fastq_wg", "fastqs_wg"), ("fastqs_wg", "fastq_wg")):
            if old in path:
                alt = path.replace(old, new)
                if os.path.isdir(alt):
                    candidates.append(alt)
    if kind == "genome":
        for old, new in (("cirulated_fasta", "circulated_fasta"), ("circulated_fasta", "cirulated_fasta")):
            if old in path:
                alt = path.replace(old, new)
                if os.path.isdir(alt):
                    candidates.append(alt)
    parent = os.path.dirname(path)
    if os.path.isdir(parent):
        base = os.path.basename(path).lower()
        try:
            for name in sorted(os.listdir(parent)):
                full = os.path.join(parent, name)
                if os.path.isdir(full) and name.lower() != base:
                    if kind == "fastq" and "fastq" in name.lower():
                        candidates.append(full)
                    elif kind == "genome" and ("fasta" in name.lower() or "genome" in name.lower()):
                        candidates.append(full)
        except OSError:
            pass
    return candidates[0] if candidates else ""


def check_existing_dir(raw: str | None, label: str = "Directory", *, kind: str = "") -> dict:
    text = str(raw or "").strip().strip('"').strip("'")
    if not text:
        return _result(False, f"Choose a {label.lower()} first.")
    path = normalize_path(text)
    if os.path.isdir(path):
        return _result(True, f"{label} OK: {path}", path=path)
    hint = suggest_dir(text, kind=kind or "fastq")
    msg = f"{label} does not exist: {path}"
    if hint:
        msg += f". Did you mean: {hint}?"
    return _result(False, msg, path=path, suggestion=hint)


def check_output_dir(raw: str | None, label: str = "Output directory") -> dict:
    """Output dirs may be created; the path or its parent must exist."""
    text = str(raw or "").strip().strip('"').strip("'")
    if not text:
        return _result(False, f"Choose a {label.lower()} first.")
    path = normalize_path(text)
    if os.path.isdir(path):
        return _result(True, f"{label} OK: {path}", path=path)
    parent = os.path.dirname(path) or os.path.dirname(path + os.sep)
    if parent and (os.path.isdir(parent) or parent in ("/", "")):
        return _result(True, f"{label} will be created at {path}", path=path, will_create=True)
    return _result(
        False,
        f"{label} cannot be created: parent folder does not exist for {path}",
        path=path,
    )


def _list_named_files(path: str, exts: tuple[str, ...], *, recurse: bool = False) -> list[str]:
    found: list[str] = []
    if not os.path.isdir(path):
        return found
    try:
        names = sorted(os.listdir(path))
    except OSError:
        return found
    for name in names:
        full = os.path.join(path, name)
        low = name.lower()
        if os.path.isfile(full) and any(low.endswith(ext) for ext in exts):
            found.append(full)
        elif recurse and os.path.isdir(full) and name not in (".", "..", ".nextflow-work", "work"):
            found.extend(_list_named_files(full, exts, recurse=False))
    return found


def check_db_dir(raw: str | None, label: str) -> dict:
    check = check_existing_dir(raw, label, kind="db")
    if not check.get("ok"):
        return check
    path = check["path"]
    try:
        entries = [n for n in os.listdir(path) if not n.startswith(".")]
    except OSError as exc:
        return _result(False, f"{label} is not readable: {path} ({exc})", path=path)
    if not entries:
        return _result(False, f"{label} is empty: {path}", path=path)
    return _result(True, f"{label} OK: {path} ({len(entries)} entr{'y' if len(entries) == 1 else 'ies'})", path=path)


def check_genome_dir(raw: str | None, extension: str = "fasta") -> dict:
    check = check_existing_dir(raw, "Genome directory", kind="genome")
    if not check.get("ok"):
        return check
    path = check["path"]
    files = _list_named_files(path, FASTA_EXTS, recurse=False)
    if not files:
        files = _list_named_files(path, FASTA_EXTS, recurse=True)
    if not files:
        ext = (extension or "fasta").lstrip(".")
        return _result(
            False,
            f"No FASTA files (*.fasta / *.fa / *.fna) found in {path}. "
            f"Expected genomes with extension like .{ext}.",
            path=path,
            files=[],
        )
    return _result(
        True,
        f"Found {len(files)} genome file(s) in {path}.",
        path=path,
        files=[os.path.basename(p) for p in files[:20]],
        file_count=len(files),
    )


def check_fastq_dir(
    raw: str | None,
    *,
    mode: str = "flye",
    concat_reads: str = "false",
    extension: str = "auto",
    short_read_dir: str = "",
) -> dict:
    """Strict FASTQ preflight (same gauge as Assembly)."""
    text = str(raw or "").strip().strip('"').strip("'")
    if not text:
        return _result(False, "Choose a FASTQ directory first.")
    path = normalize_path(text)
    if not os.path.isdir(path):
        hint = suggest_dir(text, kind="fastq")
        msg = f"FASTQ directory does not exist: {path}"
        if hint:
            msg += f". Did you mean: {hint}?"
        return _result(False, msg, path=path, suggestion=hint, files=[], concat_reads=concat_reads, extension=extension or "auto")

    try:
        from pool_reads import collect_sample_groups, inspect_layout, list_read_files
    except ImportError:
        names = [
            n
            for n in sorted(os.listdir(path))
            if n.lower().endswith(FASTQ_EXTS)
        ]
        if not names:
            return _result(False, f"No FASTQ files found in {path}.", path=path, files=[])
        return _result(
            True,
            f"Found {len(names)} FASTQ file(s) in {path}.",
            path=path,
            files=names[:20],
            concat_reads=concat_reads,
            extension=extension or "auto",
        )

    layout = inspect_layout(path, extension or "auto")
    sniffed = layout.get("extension") or extension or "auto"
    want_concat = _truthy(concat_reads) or bool(layout.get("needs_concat"))
    files: list[str] = []
    if want_concat:
        groups = collect_sample_groups(path, sniffed if sniffed != "auto" else "auto")
        files = [name for name, _srcs in groups]
        concat_reads = "true"
        if sniffed and sniffed != "auto":
            extension = sniffed
        if not files:
            return _result(
                False,
                f"No FASTQ files found under {path} for pooling. "
                "Expected sample subfolders with reads, or top-level FASTQ files.",
                path=path,
                files=[],
                concat_reads=concat_reads,
                extension=extension or sniffed,
            )
    else:
        top = layout.get("top_level_files") or []
        pooled = layout.get("pooled_files") or []
        files = top or pooled
        if not files:
            files = [os.path.basename(p) for p in list_read_files(path, recurse=False)]
        concat_reads = "false"
        if sniffed and sniffed != "auto":
            extension = sniffed
        if not files:
            return _result(
                False,
                f"No FASTQ files found in {path}. "
                "Expected files like *.fastq / *.fastq.gz "
                "(or turn Concat Reads on for barcode/sample subfolders).",
                path=path,
                files=[],
                concat_reads=concat_reads,
                extension=extension or sniffed,
            )

    if mode == "spades":
        all_names = [os.path.basename(p) for p in list_read_files(path, recurse=True)]
        lower = [n.lower() for n in all_names]
        has_pair = any("_r1" in n or n.endswith("1.fastq.gz") or n.endswith("1.fastq") for n in lower)
        if not has_pair and len(all_names) < 2:
            return _result(
                False,
                f"No Illumina paired-end FASTQs found in {path}. "
                "Expected names like sample_R1.fastq.gz / sample_R2.fastq.gz.",
                path=path,
                files=all_names[:20],
                concat_reads="false",
                extension=extension or sniffed,
            )

    if mode == "unicycler" and short_read_dir:
        short = check_fastq_dir(
            short_read_dir, mode="spades", concat_reads="false", extension="auto"
        )
        if not short.get("ok"):
            return _result(
                False,
                f"Short-read directory: {short.get('message')}",
                path=path,
                files=files[:20],
                concat_reads=concat_reads,
                extension=extension or sniffed,
            )

    return _result(
        True,
        f"Found {len(files)} sample/file(s) in {path}.",
        path=path,
        files=files[:20],
        file_count=len(files),
        concat_reads=concat_reads,
        extension=extension or sniffed,
        needs_concat=bool(layout.get("needs_concat")),
    )


def _first_failure(checks: list[dict]) -> dict | None:
    for item in checks:
        if not item.get("ok"):
            return item
    return None


def validate_form_paths(module: str, form) -> dict:
    """Validate directory fields for qc | assembly | post from a form mapping."""
    get = form.get if hasattr(form, "get") else lambda k, d="": form.get(k, d) if isinstance(form, dict) else d
    module = (module or "").strip().lower()
    fields: list[dict] = []

    def add(name: str, check: dict, *, required: bool = True):
        item = {
            "name": name,
            "ok": bool(check.get("ok")),
            "message": check.get("message") or "",
            "path": check.get("path") or check.get("directory") or "",
            "required": required,
        }
        if check.get("suggestion"):
            item["suggestion"] = check["suggestion"]
        fields.append(item)
        return check

    if module == "qc":
        add("fastq_dir", check_fastq_dir(
            get("fastq_dir", ""),
            mode="spades" if str(get("read_platform", "")).lower() == "illumina" else "flye",
            concat_reads=get("concat_reads", "false"),
            extension=(get("extension") or "auto").strip() or "auto",
        ))
        add("out_dir", check_output_dir(get("out_dir", "") or "./bactflow_out", "Output directory"))

    elif module == "assembly":
        run_flye = _truthy(get("run_flye", "false"))
        run_unicycler = _truthy(get("run_unicycler", "false"))
        run_spades = _truthy(get("run_spades", "false"))
        run_pacbio = _truthy(get("run_pacbio", "false"))
        assembler = str(get("assembler", "")).strip().lower()
        if assembler == "flye":
            run_flye = True
        elif assembler == "pacbio":
            run_pacbio = True
        elif assembler == "unicycler":
            run_unicycler = True
        elif assembler == "spades":
            run_spades = True

        needs_fastq = run_flye or run_unicycler or run_spades or run_pacbio
        if needs_fastq:
            if run_pacbio:
                mode = "pacbio"
            elif run_spades:
                mode = "spades"
            elif run_unicycler:
                mode = "unicycler"
            else:
                mode = "flye"
            add(
                "fastq_dir",
                check_fastq_dir(
                    get("fastq_dir", ""),
                    mode=mode,
                    concat_reads=get("concat_reads", "false"),
                    extension=(get("extension") or "auto").strip() or "auto",
                    short_read_dir=get("short_read_dir", "") if run_unicycler else "",
                ),
            )
            if run_unicycler:
                add(
                    "short_read_dir",
                    check_fastq_dir(
                        get("short_read_dir", ""),
                        mode="spades",
                        concat_reads="false",
                        extension="auto",
                    ),
                )
        add("out_dir", check_output_dir(get("out_dir", "") or "./bactflow_out", "Output directory"))

        if _truthy(get("bakta_annot", "false")):
            add("bakta_data_path", check_db_dir(get("bakta_data_path", ""), "Bakta database directory"))
        if _truthy(get("tax_class", "false")):
            add("gtdbtk_data_path", check_db_dir(get("gtdbtk_data_path", ""), "GTDB-Tk database directory"))
        if _truthy(get("run_checkm", "false")):
            add("checkm_db", check_db_dir(get("checkm_db", ""), "CheckM database directory"))
        if _truthy(get("run_plasmids", "false")):
            add("genomad_db", check_db_dir(get("genomad_db", ""), "geNomad database directory"))

    elif module == "post":
        setup_only = _truthy(get("setup_only", "false"))
        if not setup_only:
            add(
                "genome_dir",
                check_genome_dir(get("genome_dir", ""), get("genome_extension", "fasta") or "fasta"),
            )
        add("out_dir", check_output_dir(get("out_dir", "") or "./bactflow_out", "Output directory"))
        if _truthy(get("bakta_annot", "false")):
            add("bakta_data_path", check_db_dir(get("bakta_data_path", ""), "Bakta database directory"))
        if _truthy(get("tax_class", "false")):
            add("gtdbtk_data_path", check_db_dir(get("gtdbtk_data_path", ""), "GTDB-Tk database directory"))
        if _truthy(get("run_checkm", "false")):
            add("checkm_db", check_db_dir(get("checkm_db", ""), "CheckM database directory"))
        if _truthy(get("run_plasmids", "false")):
            add("genomad_db", check_db_dir(get("genomad_db", ""), "geNomad database directory"))
    else:
        return _result(False, f"Unknown module for path validation: {module}")

    failed = [f for f in fields if not f.get("ok")]
    if failed:
        msg = failed[0]["message"]
        if len(failed) > 1:
            msg += f" (+{len(failed) - 1} more path issue(s))"
        return {
            "ok": False,
            "message": msg,
            "fields": fields,
        }
    return {
        "ok": True,
        "message": "All directory paths look valid.",
        "fields": fields,
    }
