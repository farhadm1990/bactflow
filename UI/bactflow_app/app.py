#!/usr/bin/env python3
"""Unified BactFlow UI: welcome hub plus QC, Assembly, and Post-assembly."""

from __future__ import annotations

import importlib.util
import os
import sys
from threading import Timer

from flask import Flask, jsonify, redirect, render_template, send_from_directory
from werkzeug.middleware.dispatcher import DispatcherMiddleware

base_dir = os.path.abspath(os.path.dirname(__file__))
ui_root = os.path.abspath(os.path.join(base_dir, ".."))
if base_dir not in sys.path:
    sys.path.insert(0, base_dir)

from run_lock import snapshot  # noqa: E402


def _load_module(name: str, path: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


qc_mod = _load_module(
    "bactflow_qc_app", os.path.join(ui_root, "pre_assem_app", "pre_assembly.py")
)
assem_mod = _load_module(
    "bactflow_assembly_app", os.path.join(ui_root, "assem_app", "assembly.py")
)
post_mod = _load_module(
    "bactflow_post_app", os.path.join(ui_root, "post_assem_app", "post_assembly.py")
)

hub = Flask(
    __name__,
    template_folder=os.path.join(base_dir, "templates"),
    static_folder=os.path.join(base_dir, "static"),
)


@hub.route("/")
def home():
    return render_template("home.html")


@hub.route("/wiki")
def wiki():
    return redirect("/#wiki")


@hub.route("/run-lock")
def run_lock_status():
    return jsonify(snapshot())


@hub.route("/host-info")
def host_info():
    return jsonify(snapshot())


@hub.route("/hub/static/<path:filename>")
def hub_static(filename):
    return send_from_directory(os.path.join(base_dir, "static"), filename)


hub.wsgi_app = DispatcherMiddleware(
    hub.wsgi_app,
    {
        "/qc": qc_mod.app.wsgi_app,
        "/assembly": assem_mod.app.wsgi_app,
        "/post": post_mod.app.wsgi_app,
    },
)


def _ui_port(default=5000):
    try:
        port = int(os.environ.get("BACTFLOW_HOST_PORT") or default)
    except (TypeError, ValueError):
        port = default
    return port if 1 <= port <= 65535 else default


def _access_url(port: int | None = None) -> str:
    port = port if port is not None else _ui_port(5000)
    explicit = (os.environ.get("BACTFLOW_PUBLIC_URL") or "").strip()
    if explicit:
        return explicit if explicit.endswith("/") else explicit + "/"
    return f"http://127.0.0.1:{port}/"


def open_browser():
    try:
        from bactflow_open_browser import open_bactflow_browser
    except ImportError:
        sys.path.insert(0, os.path.join(ui_root, "assem_app"))
        from bactflow_open_browser import open_bactflow_browser
    open_bactflow_browser(port=_ui_port(5000))


if __name__ == "__main__":
    port = _ui_port(5000)
    url = _access_url(port)
    print(f"\n=== BactFlow UI ===\nOpen this URL: {url}\n", flush=True)
    want_browser = (
        os.environ.get("BACTFLOW_NO_BROWSER") != "1"
        or os.environ.get("BACTFLOW_FORCE_BROWSER") == "1"
    )
    if want_browser:
        Timer(1, open_browser).start()
    flask_debug = os.environ.get("BACTFLOW_FLASK_DEBUG") == "1"
    hub.run(
        debug=flask_debug,
        port=port,
        host="0.0.0.0",
        use_reloader=False,
        threaded=True,
    )
