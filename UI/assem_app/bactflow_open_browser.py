"""Open BactFlow UI in a browser — works for local runs and Docker on the host."""

from __future__ import annotations

import os
import socket
import subprocess
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


def _default_port() -> int:
    for key in ("BACTFLOW_HOST_PORT", "PORT"):
        raw = os.environ.get(key)
        if raw and str(raw).isdigit():
            return int(raw)
    # Infer from common module entrypoints / env
    mod = (os.environ.get("BACTFLOW_MODULE") or "").lower()
    if mod in ("preassem", "pre", "pre_assembly"):
        return 5000
    if mod in ("assem", "assembly"):
        return 5002
    if mod in ("postassem", "post", "post_assembly"):
        return 5001
    return 5000


def public_ui_url(port: int | None = None) -> str:
    explicit = (os.environ.get("BACTFLOW_PUBLIC_URL") or "").strip()
    if explicit:
        return explicit if explicit.endswith("/") else explicit + "/"
    host = (os.environ.get("BACTFLOW_PUBLIC_HOST") or "127.0.0.1").strip() or "127.0.0.1"
    p = port or _default_port()
    return f"http://{host}:{p}/"


def _in_docker() -> bool:
    return os.environ.get("BACTFLOW_IN_DOCKER") == "1" or Path("/.dockerenv").exists()


def _write_url_dropfile(url: str) -> None:
    """Write URL where a host-side hook (watching ~/.bactflow/open_url) can see it."""
    candidates = []
    # Raw docker often mounts host $HOME to /home
    for base in (
        Path("/home/.bactflow"),
        Path("/home/farhad/.bactflow"),
        Path(os.path.expanduser("~/.bactflow")),
    ):
        candidates.append(base)
    home = os.environ.get("HOME")
    if home:
        candidates.append(Path(home) / ".bactflow")

    seen = set()
    for folder in candidates:
        key = str(folder)
        if key in seen:
            continue
        seen.add(key)
        try:
            folder.mkdir(parents=True, exist_ok=True)
            target = folder / "open_url"
            target.write_text(url + "\n", encoding="utf-8")
        except OSError:
            continue


def _hook_targets(port: int = 19264) -> list[str]:
    hosts = []
    env_host = (os.environ.get("BACTFLOW_BROWSER_HOOK_HOST") or "").strip()
    if env_host:
        hosts.append(env_host)
    # Prefer docker bridge gateway first — host.docker.internal often needs --add-host.
    hosts.append("172.17.0.1")
    try:
        with open("/proc/net/route", "r", encoding="utf-8") as handle:
            next(handle, None)
            for line in handle:
                parts = line.split()
                if len(parts) >= 3 and parts[1] == "00000000":
                    raw = parts[2]
                    # little-endian hex gateway
                    gw = socket.inet_ntoa(bytes.fromhex(raw)[::-1])
                    hosts.append(gw)
                    break
    except OSError:
        pass
    hosts.append("host.docker.internal")
    # de-dupe
    out = []
    seen = set()
    for h in hosts:
        if h and h not in seen:
            seen.add(h)
            out.append(f"http://{h}:{port}/open")
    return out


def _request_host_hook(url: str) -> bool:
    port = int(os.environ.get("BACTFLOW_BROWSER_HOOK_PORT") or 19264)
    payload = urllib.parse.urlencode({"url": url}).encode("utf-8")
    for endpoint in _hook_targets(port):
        try:
            req = urllib.request.Request(
                endpoint + "?" + urllib.parse.urlencode({"url": url}),
                data=payload,
                method="POST",
                headers={"Content-Type": "application/x-www-form-urlencoded"},
            )
            with urllib.request.urlopen(req, timeout=1.5) as resp:
                if 200 <= getattr(resp, "status", 200) < 300:
                    return True
        except (urllib.error.URLError, TimeoutError, OSError):
            try:
                with urllib.request.urlopen(
                    endpoint + "?" + urllib.parse.urlencode({"url": url}),
                    timeout=1.5,
                ) as resp:
                    if 200 <= getattr(resp, "status", 200) < 300:
                        return True
            except (urllib.error.URLError, TimeoutError, OSError):
                continue
    return False


def _local_xdg_open(url: str) -> bool:
    env = os.environ.copy()
    if not env.get("DISPLAY"):
        if Path("/tmp/.X11-unix/X0").exists():
            env["DISPLAY"] = ":0"
        elif Path("/tmp/.X11-unix/X1").exists():
            env["DISPLAY"] = ":1"
    for cmd in (
        ["xdg-open", url],
        ["gio", "open", url],
        ["firefox", url],
        ["google-chrome", url],
        ["chromium-browser", url],
    ):
        try:
            subprocess.Popen(
                cmd,
                env=env,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                start_new_session=True,
            )
            return True
        except OSError:
            continue
    return False


def open_bactflow_browser(port: int | None = None, delay: float = 0.0) -> str:
    """Open the UI URL. Returns the URL attempted."""
    if delay > 0:
        time.sleep(delay)
    url = public_ui_url(port)

    if os.environ.get("BACTFLOW_NO_BROWSER") == "1" and os.environ.get("BACTFLOW_FORCE_BROWSER") != "1":
        # Host launcher (bactflow.sh) owns the browser open in this mode.
        return url

    # 1) Host hook (best for Docker / raw docker run)
    if _in_docker():
        _write_url_dropfile(url)
        if _request_host_hook(url):
            return url

    # 2) Direct local open (native runs or X11-forwarded containers)
    try:
        import webbrowser

        webbrowser.open(url, new=2, autoraise=True)
    except Exception:
        pass
    _local_xdg_open(url)
    return url
