"""Parse and enforce a per-run RAM cap so heavy tools cannot OOM the host.

Nextflow process.memory alone is only a scheduler hint. This module also:
  - respects BACTFLOW_MEMORY / BACTFLOW_DOCKER_MEMORY from bactflow.sh --memory
  - wraps the pipeline shell in a cgroup (systemd-run) when available
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess


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


def _parse_size_to_gb(raw) -> int | None:
    text = str(raw or "").strip().lower().replace(" ", "")
    if not text or text.startswith("auto"):
        return None
    text = text.replace("gi", "g").replace("gb", "g")
    mult = 1.0
    if text.endswith("t"):
        mult = 1024.0
        text = text[:-1]
    elif text.endswith("g"):
        mult = 1.0
        text = text[:-1]
    elif text.endswith("m"):
        mult = 1.0 / 1024.0
        text = text[:-1]
    elif text.endswith("k"):
        mult = 1.0 / (1024.0 * 1024.0)
        text = text[:-1]
    try:
        return max(1, int(round(float(text) * mult)))
    except (TypeError, ValueError):
        return None


def env_launch_cap_gb() -> int | None:
    """Hard ceiling from bactflow.sh --memory (exported into the UI process)."""
    for key in (
        "BACTFLOW_MEMORY_GB",
        "BACTFLOW_DOCKER_MEMORY",
        "BACTFLOW_MEMORY",
    ):
        gb = _parse_size_to_gb(os.environ.get(key))
        if gb:
            return gb
    return None


def default_memory_gb(host: int | None = None) -> int:
    host = int(host or host_mem_gb())
    launch = env_launch_cap_gb()
    reserve = 1 if host > 3 else 0
    base = max(2, min(8, host - reserve))
    if launch:
        return max(2, min(base, launch))
    return base


def parse_and_cap(raw, host: int | None = None):
    """Return (gb, was_capped, host_gb, max_allowed)."""
    host = int(host or host_mem_gb())
    reserve = 1 if host > 3 else 0
    max_allowed = max(2, host - reserve)
    launch = env_launch_cap_gb()
    if launch:
        max_allowed = min(max_allowed, launch)

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


_SYSTEMD_OK: bool | None = None


def _systemd_run_works() -> bool:
    global _SYSTEMD_OK
    if _SYSTEMD_OK is not None:
        return _SYSTEMD_OK
    if not shutil.which("systemd-run"):
        _SYSTEMD_OK = False
        return False
    try:
        # Probe a tiny scope; fails on hosts without a user systemd session.
        proc = subprocess.run(
            [
                "systemd-run",
                "--user",
                "--quiet",
                "--collect",
                "--scope",
                "-p",
                "MemoryMax=64M",
                "-p",
                "MemorySwapMax=0",
                "--",
                "true",
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=8,
            check=False,
        )
        _SYSTEMD_OK = proc.returncode == 0
    except (OSError, subprocess.TimeoutExpired):
        _SYSTEMD_OK = False
    return _SYSTEMD_OK


def hard_cap_wrapper(mem_gb: int) -> tuple[str, str]:
    """Return (shell_prefix, mode) that wraps the rest of a bash -c script.

    mode is 'systemd-user', 'docker-cgroup', or 'soft-only'.
    """
    mem_gb = max(2, int(mem_gb))
    # Inside Docker / Apptainer the container cgroup already hard-caps when
    # bactflow.sh passed --memory (+ --memory-swap).
    if os.environ.get("BACTFLOW_IN_DOCKER") == "1" or os.path.exists("/.dockerenv"):
        return "", "docker-cgroup"
    if _systemd_run_works():
        # Wrap the remainder of the script under MemoryMax. The pipeline is
        # invoked as: bash -c '<preamble><wrapper> bash -c <quoted pipeline>'
        return (
            f'systemd-run --user --quiet --collect --scope '
            f'-p MemoryMax={mem_gb}G -p MemorySwapMax=0 -- bash -c ',
            "systemd-user",
        )
    return "", "soft-only"


def inject_nextflow_memory(command: str, out_dir: str, cpus, mem_gb: int) -> str:
    cfg = write_resource_config(out_dir, cpus, mem_gb)
    extra = f"-c '{cfg}' --memory_gb {int(mem_gb)} "
    if "-ansi-log" in command:
        command = command.replace("-ansi-log", extra + "-ansi-log", 1)
    else:
        command = command.rstrip() + " " + extra

    prefix, mode = hard_cap_wrapper(mem_gb)
    if prefix and mode == "systemd-user":
        # Quote the whole pipeline for the inner bash -c.
        import shlex

        command = prefix + shlex.quote(command)
    return command


_MEM_NOTE_RE = re.compile(r"Memory cap:", re.I)


def preamble(mem_gb: int, host: int, capped: bool) -> str:
    launch = env_launch_cap_gb()
    _, mode = hard_cap_wrapper(mem_gb)
    if mode == "docker-cgroup":
        enforce = "container cgroup (hard)"
    elif mode == "systemd-user":
        enforce = "systemd MemoryMax (hard)"
    else:
        enforce = "Nextflow only (soft — start with docker --memory for a hard OS cap)"

    note = (
        f"echo 'Memory cap: {int(mem_gb)} GB on a {int(host)} GB host "
        f"[{enforce}].'"
    )
    if launch:
        note += f"\necho 'Launcher --memory ceiling: {int(launch)} GB.'"
    if capped:
        note += (
            f"\necho 'Requested RAM was lowered to {int(mem_gb)} GB "
            "so this machine / launcher limit is not exceeded.'"
        )
    if mode == "soft-only":
        note += (
            "\necho 'WARNING: no hard cgroup available; RAM may exceed the "
            "form value. Prefer: ./bactflow.sh run … --service docker --memory "
            f"{int(mem_gb)}g'"
        )
    return note + "\n"
