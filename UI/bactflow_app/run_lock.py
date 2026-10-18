"""Cross-module run lock so only one BactFlow process runs at a time."""

from __future__ import annotations

import json
import os
import signal
import time

LOCK_DIR = os.path.join(os.environ.get("HOME") or "/tmp", ".bactflow")
LOCK_PATH = os.path.join(LOCK_DIR, "ui-run.lock")

# Sync locks use the Flask PID; if a request dies without release they can stick.
# Drop them after this many seconds so the hub is not blocked forever.
SYNC_MAX_AGE_SEC = 45 * 60

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
    # Always drop locks whose owner PID is gone (previous bug: only "worker" was checked,
    # so a dead sync lock after a restart kept the UI stuck on "QC is running").
    if not _pid_alive(data.get("pid")):
        _clear()
        return None
    if kind == "sync":
        started = data.get("updated") or data.get("started") or 0
        try:
            age = time.time() - float(started)
        except (TypeError, ValueError):
            age = SYNC_MAX_AGE_SEC + 1
        if age > SYNC_MAX_AGE_SEC:
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


def _terminate_tree(pid: int) -> list[int]:
    """Best-effort terminate of pid and children. Never targets the current process."""
    killed: list[int] = []
    if pid <= 0 or pid == os.getpid():
        return killed
    try:
        import psutil

        parent = psutil.Process(pid)
        procs = list(parent.children(recursive=True)) + [parent]
        for proc in procs:
            try:
                if proc.pid == os.getpid():
                    continue
                proc.terminate()
                killed.append(proc.pid)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue
        _, alive = psutil.wait_procs(procs, timeout=2)
        for proc in alive:
            try:
                if proc.pid != os.getpid():
                    proc.kill()
                    if proc.pid not in killed:
                        killed.append(proc.pid)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue
        return killed
    except Exception:
        pass
    try:
        os.kill(pid, signal.SIGTERM)
        killed.append(pid)
    except (OSError, ProcessLookupError):
        pass
    return killed


def force_reset() -> dict:
    """Clear the run lock and stop any worker PID it pointed at (not the UI itself)."""
    data = _read()
    killed: list[int] = []
    module = (data or {}).get("module")
    label = (data or {}).get("label")
    pid = (data or {}).get("pid")
    kind = (data or {}).get("kind")
    try:
        pid_i = int(pid) if pid is not None else 0
    except (TypeError, ValueError):
        pid_i = 0
    if data and pid_i > 0 and pid_i != os.getpid() and _pid_alive(pid_i):
        killed = _terminate_tree(pid_i)
    _clear()
    return {
        "had_lock": bool(data),
        "module": module,
        "label": label,
        "pid": pid_i or None,
        "kind": kind,
        "killed": killed,
        "cleared": True,
        "lock_path": LOCK_PATH,
    }


def snapshot() -> dict:
    from memory_cap import default_memory_gb, env_launch_cap_gb, host_mem_gb

    host = host_mem_gb()
    launch = env_launch_cap_gb()
    max_gb = max(2, host - (1 if host > 3 else 0))
    if launch:
        max_gb = min(max_gb, launch)
    data = current()
    return {
        "running": bool(data),
        "module": (data or {}).get("module"),
        "label": (data or {}).get("label"),
        "pid": (data or {}).get("pid"),
        "kind": (data or {}).get("kind"),
        "lock_path": LOCK_PATH,
        "host_mem_gb": host,
        "launch_memory_gb": launch,
        "max_memory_gb": max_gb,
        "default_memory_gb": default_memory_gb(host),
    }
