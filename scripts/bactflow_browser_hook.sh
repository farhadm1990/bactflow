#!/usr/bin/env bash
# Host-side helper: open BactFlow UI URLs in the desktop browser.
# Used when modules run in Docker (raw `docker run` or ./bactflow.sh).
#
# Usage:
#   bactflow_browser_hook.sh --install   # start daemon + autostart
#   bactflow_browser_hook.sh --daemon    # run in foreground
#   bactflow_browser_hook.sh --once URL  # open one URL and exit
#   bactflow_browser_hook.sh --status

set -euo pipefail

HOOK_PORT="${BACTFLOW_BROWSER_HOOK_PORT:-19264}"
STATE_DIR="${BACTFLOW_HOOK_STATE_DIR:-${HOME}/.bactflow}"
URL_FILE="${STATE_DIR}/open_url"
PID_FILE="${STATE_DIR}/browser_hook.pid"
LOG_FILE="${STATE_DIR}/browser_hook.log"
AUTOSTART="${HOME}/.config/autostart/bactflow-browser-hook.desktop"

mkdir -p "${STATE_DIR}"

open_url() {
    local url="$1"
    [[ -n "${url}" ]] || return 1
    echo "[bactflow-hook] Opening ${url}" >>"${LOG_FILE}"

    if [[ -z "${DISPLAY:-}" ]]; then
        if [[ -S /tmp/.X11-unix/X0 ]]; then
            export DISPLAY=:0
        elif [[ -S /tmp/.X11-unix/X1 ]]; then
            export DISPLAY=:1
        fi
    fi
    if [[ -z "${XDG_RUNTIME_DIR:-}" && -d "/run/user/$(id -u)" ]]; then
        export XDG_RUNTIME_DIR="/run/user/$(id -u)"
    fi
    if [[ -z "${DBUS_SESSION_BUS_ADDRESS:-}" && -n "${XDG_RUNTIME_DIR:-}" && -S "${XDG_RUNTIME_DIR}/bus" ]]; then
        export DBUS_SESSION_BUS_ADDRESS="unix:path=${XDG_RUNTIME_DIR}/bus"
    fi

    if command -v python3 >/dev/null 2>&1; then
        python3 - "${url}" <<'PY' >>"${LOG_FILE}" 2>&1 || true
import sys, webbrowser
webbrowser.open(sys.argv[1], new=2, autoraise=True)
PY
        return 0
    fi

    local opener
    for opener in xdg-open gio open gnome-open sensible-browser firefox google-chrome chromium-browser chromium; do
        if command -v "${opener}" >/dev/null 2>&1; then
            if [[ "${opener}" == "gio" ]]; then
                nohup gio open "${url}" >>"${LOG_FILE}" 2>&1 &
            else
                nohup "${opener}" "${url}" >>"${LOG_FILE}" 2>&1 &
            fi
            return 0
        fi
    done
    echo "[bactflow-hook] No browser opener found for ${url}" >>"${LOG_FILE}"
    return 1
}

is_running() {
    if [[ -f "${PID_FILE}" ]]; then
        local pid
        pid="$(cat "${PID_FILE}" 2>/dev/null || true)"
        if [[ -n "${pid}" ]] && kill -0 "${pid}" 2>/dev/null; then
            return 0
        fi
    fi
    return 1
}

run_daemon() {
    if is_running; then
        echo "BactFlow browser hook already running (pid $(cat "${PID_FILE}"))"
        return 0
    fi

    echo $$ >"${PID_FILE}"
    echo "[bactflow-hook] Listening on 0.0.0.0:${HOOK_PORT} and watching ${URL_FILE}" >>"${LOG_FILE}"

    # Prefer Python HTTP server (works without extra packages).
    exec python3 - "${HOOK_PORT}" "${URL_FILE}" "${LOG_FILE}" <<'PY'
import os, sys, time, threading, webbrowser
from http.server import BaseHTTPRequestHandler, HTTPServer
from urllib.parse import urlparse, parse_qs, unquote_plus

port = int(sys.argv[1])
url_file = sys.argv[2]
log_file = sys.argv[3]

def log(msg):
    with open(log_file, "a", encoding="utf-8") as fh:
        fh.write(msg.rstrip() + "\n")

def open_url(url: str):
    if not url:
        return
    log(f"[bactflow-hook] Opening {url}")
    try:
        webbrowser.open(url, new=2, autoraise=True)
    except Exception as exc:
        log(f"[bactflow-hook] webbrowser failed: {exc}")
        for cmd in ("xdg-open", "gio", "firefox", "google-chrome", "chromium-browser"):
            if cmd == "gio":
                os.spawnlp(os.P_NOWAIT, "gio", "gio", "open", url)
            else:
                try:
                    os.spawnlp(os.P_NOWAIT, cmd, cmd, url)
                    break
                except Exception:
                    continue

class Handler(BaseHTTPRequestHandler):
    def log_message(self, fmt, *args):
        log("[bactflow-hook-http] " + (fmt % args))

    def _send(self, code=200, body=b"ok"):
        self.send_response(code)
        self.send_header("Content-Type", "text/plain")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        parsed = urlparse(self.path)
        if parsed.path in ("/health", "/"):
            self._send(200, b"ok")
            return
        if parsed.path == "/open":
            qs = parse_qs(parsed.query)
            url = (qs.get("url") or [""])[0]
            open_url(unquote_plus(url))
            self._send(200, b"opened")
            return
        self._send(404, b"not found")

    def do_POST(self):
        length = int(self.headers.get("Content-Length") or 0)
        raw = self.rfile.read(length).decode("utf-8", errors="replace") if length else ""
        parsed = urlparse(self.path)
        qs = parse_qs(parsed.query)
        url = (qs.get("url") or [""])[0]
        if not url and raw:
            # url=... or plain body
            if raw.startswith("url="):
                url = unquote_plus(raw.split("=", 1)[1].strip())
            else:
                url = raw.strip()
        if parsed.path == "/open" and url:
            open_url(url)
            self._send(200, b"opened")
            return
        self._send(400, b"missing url")

def watch_file():
    last = None
    while True:
        try:
            if os.path.isfile(url_file):
                with open(url_file, "r", encoding="utf-8") as fh:
                    url = fh.read().strip()
                if url and url != last:
                    open_url(url)
                    last = url
                    try:
                        os.remove(url_file)
                    except OSError:
                        pass
        except Exception as exc:
            log(f"[bactflow-hook] watch error: {exc}")
        time.sleep(0.7)

threading.Thread(target=watch_file, daemon=True).start()
httpd = HTTPServer(("0.0.0.0", port), Handler)
log(f"[bactflow-hook] HTTP server on 0.0.0.0:{port}")
httpd.serve_forever()
PY
}

install_hook() {
    mkdir -p "$(dirname "${AUTOSTART}")"
    local self
    self="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
    cat >"${AUTOSTART}" <<EOF
[Desktop Entry]
Type=Application
Name=BactFlow Browser Hook
Comment=Open BactFlow Docker UI URLs in your browser
Exec=${self} --daemon
X-GNOME-Autostart-enabled=true
StartupNotify=false
Terminal=false
EOF

    if is_running; then
        echo "Browser hook already running."
    else
        nohup "$self" --daemon >>"${LOG_FILE}" 2>&1 &
        sleep 0.4
        if is_running; then
            echo "Started BactFlow browser hook on port ${HOOK_PORT}."
        else
            echo "Warning: failed to confirm hook start; check ${LOG_FILE}" >&2
        fi
    fi
    echo "Autostart: ${AUTOSTART}"
}

case "${1:-}" in
    --install)
        install_hook
        ;;
    --daemon)
        run_daemon
        ;;
    --once)
        open_url "${2:-}"
        ;;
    --status)
        if is_running; then
            echo "running pid=$(cat "${PID_FILE}") port=${HOOK_PORT}"
        else
            echo "stopped"
            exit 1
        fi
        ;;
    --stop)
        if [[ -f "${PID_FILE}" ]]; then
            kill "$(cat "${PID_FILE}")" 2>/dev/null || true
            rm -f "${PID_FILE}"
        fi
        echo "stopped"
        ;;
    *)
        cat <<EOF
Usage: $0 --install | --daemon | --once URL | --status | --stop
EOF
        exit 1
        ;;
esac
