"""Cross-module run lock so only one BactFlow process runs at a time."""

from __future__ import annotations

import json
import os
import time

LOCK_DIR = os.path.join(os.environ.get("HOME") or "/tmp", ".bactflow")
LOCK_PATH = os.path.join(LOCK_DIR, "ui-run.lock")

LABELS = {
    "qc": "QC",
    "assembly": "Assembly",
    "post": "Post-assembly",
}


def _pid_alive(pid) -> bool:
    try:
        pid = int(pid)
    except (TypeError, ValueError):
        return False
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
        return True
    except (OSError, ProcessLookupError):
        return False


def _read() -> dict | None:
    try:
        with open(LOCK_PATH, "r", encoding="utf-8") as handle:
            data = json.load(handle)
        if isinstance(data, dict) and data.get("module"):
            return data
    except (OSError, json.JSONDecodeError, TypeError):
        return None
    return None


def _write(data: dict) -> None:
    os.makedirs(LOCK_DIR, exist_ok=True)
    tmp = LOCK_PATH + ".tmp"
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(data, handle)
        handle.write("\n")
    os.replace(tmp, LOCK_PATH)


def _clear() -> None:
    try:
        os.remove(LOCK_PATH)
    except OSError:
        pass


def current() -> dict | None:
    data = _read()
    if not data:
        return None
    kind = data.get("kind") or "worker"
    if kind == "worker" and not _pid_alive(data.get("pid")):
        _clear()
        return None
    return data


def try_acquire(module: str, pid: int, label: str | None = None, kind: str = "sync"):
    """Claim the global slot. Returns (ok, error_message)."""
    existing = current()
    if existing:
        other = existing.get("module")
        other_label = existing.get("label") or LABELS.get(other, other)
        if other != module:
            return False, (
                f"{other_label} is still running. Finish or stop it before launching another BactFlow process."
            )
        return False, f"{other_label} is already running."
    payload = {
        "module": module,
        "label": label or LABELS.get(module, module),
        "pid": int(pid) if pid else os.getpid(),
        "kind": kind,
        "started": time.time(),
    }
    _write(payload)
    return True, ""


def update_lock(pid=None, kind=None, label=None):
    data = _read() or {}
    if pid is not None:
        data["pid"] = int(pid)
    if kind:
        data["kind"] = kind
    if label:
        data["label"] = label
    data["updated"] = time.time()
    if data.get("module"):
        _write(data)


def release(module: str | None = None, pid=None) -> None:
    data = _read()
    if not data:
        return
    if module and data.get("module") != module:
        return
    if pid is not None:
        try:
            if int(data.get("pid") or 0) != int(pid):
                return
        except (TypeError, ValueError):
            return
    _clear()


def snapshot() -> dict:
    from memory_cap import default_memory_gb, host_mem_gb

    host = host_mem_gb()
    data = current()
    return {
        "running": bool(data),
        "module": (data or {}).get("module"),
        "label": (data or {}).get("label"),
        "pid": (data or {}).get("pid"),
        "host_mem_gb": host,
        "max_memory_gb": max(2, host - (1 if host > 3 else 0)),
        "default_memory_gb": default_memory_gb(host),
    }
