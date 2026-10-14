"""Parse and enforce a per-run RAM cap so heavy tools cannot OOM the host."""

from __future__ import annotations

import os
import re


def host_mem_gb() -> int:
    try:
        with open("/proc/meminfo", "r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("MemTotal:"):
                    kb = int(line.split()[1])
                    return max(1, kb // (1024 * 1024))
    except (OSError, ValueError, IndexError):
        pass
    try:
        import psutil

        return max(1, int(psutil.virtual_memory().total // (1024 ** 3)))
    except Exception:
        return 8


def default_memory_gb(host: int | None = None) -> int:
    host = int(host or host_mem_gb())
    reserve = 1 if host > 3 else 0
    return max(2, min(8, host - reserve))


def parse_and_cap(raw, host: int | None = None):
    """Return (gb, was_capped, host_gb, max_allowed)."""
    host = int(host or host_mem_gb())
    reserve = 1 if host > 3 else 0
    max_allowed = max(2, host - reserve)
    text = str(raw or "").strip().lower().replace("gb", "").replace("g", "").strip()
    try:
        gb = int(float(text))
    except (TypeError, ValueError):
        gb = default_memory_gb(host)
    gb = max(2, gb)
    capped = gb > max_allowed
    if capped:
        gb = max_allowed
    return gb, capped, host, max_allowed


def nxf_opts_export(mem_gb: int) -> str:
    heap = int(int(mem_gb) * 1024 * 0.25)
    heap = max(512, min(4096, heap))
    return (
        f'export NXF_OPTS="-Xmx{heap}m"\n'
        f"export BACTFLOW_MEMORY_GB={int(mem_gb)}\n"
    )


def write_resource_config(out_dir: str, cpus, mem_gb: int) -> str:
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, ".bactflow-resources.config")
    try:
        cpus_n = max(1, int(cpus or 1))
    except (TypeError, ValueError):
        cpus_n = 1
    mem_gb = max(2, int(mem_gb))
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(
            f"params.memory = '{mem_gb} GB'\n"
            f"params.memory_gb = {mem_gb}\n"
            "executor {\n"
            f"    cpus = {cpus_n}\n"
            f"    memory = '{mem_gb} GB'\n"
            "}\n"
            "process {\n"
            "    memory = params.memory\n"
            "}\n"
        )
    return path


def inject_nextflow_memory(command: str, out_dir: str, cpus, mem_gb: int) -> str:
    cfg = write_resource_config(out_dir, cpus, mem_gb)
    extra = f"-c '{cfg}' --memory_gb {int(mem_gb)} "
    if "-ansi-log" in command:
        return command.replace("-ansi-log", extra + "-ansi-log", 1)
    return command.rstrip() + " " + extra


_MEM_NOTE_RE = re.compile(r"Memory cap:", re.I)


def preamble(mem_gb: int, host: int, capped: bool) -> str:
    note = f"echo 'Memory cap: {int(mem_gb)} GB on a {int(host)} GB host.'"
    if capped:
        note += (
            f"\necho 'Requested RAM was lowered to {int(mem_gb)} GB "
            "so this machine is not exhausted.'"
        )
    return note + "\n"
