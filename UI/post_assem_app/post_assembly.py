#!/usr/bin/env python3
import sys, subprocess, importlib, os, shlex
import pty
import select
import hashlib
from flask import Flask, render_template, request, redirect, Response, send_from_directory, stream_with_context, jsonify, render_template_string, abort, url_for
import json
from datetime import datetime, timezone
import pandas as pd
import numpy as np
import os
import subprocess
import sys
from multiprocessing import Process, Manager, Queue
import time
import signal
import psutil  # for better process control

psutil.cpu_percent(interval=None)
_JOB_CPU_CACHE = {"pid": None, "stamp": 0.0, "cpu_sec": 0.0}
from threading import Timer
import webbrowser
import threading
import base64
import glob
import re


base_dir = os.path.abspath(os.path.dirname(__file__))# we can have access to all files from everywhere
IN_DOCKER = os.environ.get("BACTFLOW_IN_DOCKER") == "1" or os.path.exists("/.dockerenv")
app = Flask(__name__, 
            template_folder = os.path.join(base_dir, "templates"),
            static_folder = os.path.join(base_dir, "static"))

sys.path.insert(0, base_dir)
_HUB = os.path.abspath(os.path.join(base_dir, "..", "bactflow_app"))
if _HUB not in sys.path:
    sys.path.insert(0, _HUB)

BF_MODULE = "post"


def _load_collect_variant_tables():
    try:
        from vcf_viewer import collect_variant_tables as fn
        return fn
    except Exception:
        pass
    try:
        import importlib.util

        path = os.path.join(base_dir, "vcf_viewer.py")
        spec = importlib.util.spec_from_file_location("bactflow_vcf_viewer", path)
        if spec is None or spec.loader is None:
            return None
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return getattr(mod, "collect_variant_tables", None)
    except Exception:
        return None


collect_variant_tables = _load_collect_variant_tables()


BACTFLOW_RUNTIME_SH = os.path.join(base_dir, "bactflow_runtime.sh")

NF_JAVA_SETUP = f"""
source "{BACTFLOW_RUNTIME_SH}"
bactflow_prepare_nextflow || exit 1
echo "Using Java for Nextflow: $JAVA_CMD"
"$JAVA_CMD" -version
echo "Using Nextflow: $BACTFLOW_NEXTFLOW_BIN"
"$BACTFLOW_NEXTFLOW_BIN" -version
export NXF_ANSI_LOG=true
"""


def with_nextflow_java(command, mem_gb=None):
    extra = ""
    if mem_gb:
        from memory_cap import nxf_opts_export
        extra = nxf_opts_export(mem_gb)
    return extra + NF_JAVA_SETUP + "\n" + command


def _lock_acquire(label, kind="sync"):
    from run_lock import try_acquire
    return try_acquire(BF_MODULE, os.getpid(), label, kind=kind)


def _lock_worker(pid):
    from run_lock import update_lock
    update_lock(pid=pid, kind="worker")


def _lock_release(pid=None):
    from run_lock import release
    release(module=BF_MODULE, pid=pid)


def _memory_from_form():
    from memory_cap import parse_and_cap
    return parse_and_cap(request.form.get("memory_gb"))


def _truthy(value, default=False):
    if value is None or str(value).strip() == "":
        return default
    return str(value).strip().lower() in ("true", "1", "yes", "on")


def _nf_work_dir(out_dir):
    return os.path.join(out_dir, ".nextflow-work")

ANSI_RE = re.compile(r"\x1b\[[0-9;?]*[A-Za-z]|\x1b\][^\x07]*\x07")


def clean_log_line(line):
    line = ANSI_RE.sub("", line or "")
    return line.replace("\r", "").strip()


def _run_bash(script):
    return subprocess.run(
        ["bash", "-lc", script],
        check=False,
        text=True,
        capture_output=True,
    )


def _first_existing(*paths):
    for path in paths:
        if path and os.path.isfile(path):
            return path
    return None


def _image_data_uri(path):
    with open(path, "rb") as img_file:
        img_base = base64.b64encode(img_file.read()).decode("utf-8")
    lower = path.lower()
    if lower.endswith(".png"):
        mime = "image/png"
    elif lower.endswith(".jpg") or lower.endswith(".jpeg"):
        mime = "image/jpeg"
    else:
        mime = "application/octet-stream"
    return f"data:{mime};base64,{img_base}"


def _selected_gene_types():
    types = [t.strip() for t in request.form.getlist("gene_type") if t.strip()]
    if not types:
        single = (request.form.get("gene_type") or "cds").strip()
        types = [single] if single else ["cds"]
    return types





#assembly
@app.route('/',  methods =['GET', 'POST'])
def assembly():

    return render_template('post-assembly.html')


def _ui_port(default=5001):
    try:
        port = int(os.environ.get("BACTFLOW_HOST_PORT") or default)
    except (TypeError, ValueError):
        port = default
    return port if 1 <= port <= 65535 else default


def open_browser():
    from bactflow_open_browser import open_bactflow_browser
    open_bactflow_browser(port=_ui_port(5001))


########################################################
#                                                      #
#                  Running assembly                    #
#                                                      #
########################################################

# define function for bactflow run in multiprocess for constant realtime steraming.

manager = Manager()
process_status = manager.dict({"pid": None, "running": False, "exit_code": None})
output_queue = Queue()
output_history = manager.list() # to store output history
STOP_USER_MSG = "Bactflow has been stopped by the user!"


def _repair_stdio_for_fork():
    """Process.start() flushes stdout/stderr; a closed pipe raises BrokenPipeError."""
    for name in ("stdout", "stderr"):
        stream = getattr(sys, name, None)
        if stream is None:
            continue
        try:
            stream.flush()
        except (BrokenPipeError, OSError):
            try:
                setattr(
                    sys,
                    name,
                    open(os.devnull, "w", encoding="utf-8", errors="replace"),
                )
            except OSError:
                pass


def _start_worker_process(proc):
    """Start a multiprocessing worker without dying on a broken stdout pipe."""
    _repair_stdio_for_fork()
    try:
        proc.start()
    except BrokenPipeError:
        _repair_stdio_for_fork()
        proc.start()


_TASK_DONE_RE = re.compile(
    r"Task completed > TaskHandler\[.*?name: (?P<name>.+?); status: (?P<status>\w+); exit: (?P<exit>[^;]+);.*?workDir: (?P<wd>[^\]\s]+)"
)


def _task_completed_line(line):
    """Turn a DEBUG 'Task completed' record into the progress line the UI parses."""
    m = _TASK_DONE_RE.search(line or "")
    if not m or m.group("status") != "COMPLETED" or m.group("exit").strip() != "0":
        return None
    parts = m.group("wd").rstrip("/").split("/")
    if len(parts) < 2:
        return None
    task_hash = f"{parts[-2]}/{parts[-1][:6]}"
    return f"[{task_hash}] process > {m.group('name').strip()} [100%] 1 of 1 \u2714"


def _file_identity(path):
    try:
        st = os.stat(path)
    except OSError:
        return None
    return (st.st_dev, st.st_ino)


def _tail_nextflow_log(work_root, output_queue, stop_event, seen_lines):
    """Forward process status/completion lines from .nextflow.log (PTY often misses them)."""
    log_path = os.path.join(work_root, ".nextflow.log")
    # Nextflow renames the previous run's log at startup; only follow the new file.
    stale = _file_identity(log_path)
    current = None
    pos = 0
    skip_re = re.compile(
        r"\[Task monitor\]|TaskPollingMonitor|\bDEBUG\b",
        re.IGNORECASE,
    )
    # Prefer the compact Nextflow status fragment when the line is a full logger record.
    status_re = re.compile(
        r"(\[[0-9a-f/]+\]\s+(?:Submitted|Cached)\s+process\s+>\s+\S+(?:\s+\([^)]*\))?)",
        re.IGNORECASE,
    )
    done_re = re.compile(
        r"✔|Executor finished|Task completed|COMPLETED|\[100%\].*of",
        re.IGNORECASE,
    )
    while not stop_event.is_set():
        ident = _file_identity(log_path)
        if ident is not None and ident != stale:
            if ident != current:
                current = ident
                pos = 0
            try:
                with open(log_path, "rb") as handle:
                    if os.fstat(handle.fileno()).st_size < pos:
                        pos = 0
                    handle.seek(pos)
                    data = handle.read()
            except OSError:
                data = b""
            end = data.rfind(b"\n") + 1
            pos += end
            chunk = data[:end].decode("utf-8", errors="replace")
            for raw in chunk.splitlines():
                line = clean_log_line(raw)
                if not line or line in seen_lines:
                    continue
                done_line = _task_completed_line(line)
                if done_line:
                    if done_line not in seen_lines:
                        seen_lines.add(line)
                        seen_lines.add(done_line)
                        output_queue.put(done_line)
                    continue
                if skip_re.search(line):
                    continue

                forward = None
                status_m = status_re.search(line)
                if status_m:
                    forward = status_m.group(1)
                elif done_re.search(line) or (
                    re.search(r"\[[0-9a-f/]+\]", line, re.I)
                    and re.search(r"process\s+>", line, re.I)
                    and re.search(r"\[\s*\d+%\]", line)
                ):
                    # Strip logger prefix when present: "... INFO ... - [ab/cd] process > ..."
                    hash_m = re.search(
                        r"(\[[0-9a-f/]+\]\s+.*)$",
                        line,
                        re.I,
                    )
                    forward = hash_m.group(1) if hash_m else line

                if forward and forward not in seen_lines:
                    seen_lines.add(line)
                    seen_lines.add(forward)
                    output_queue.put(forward)
        stop_event.wait(0.8)


def _pipeline_alive(pid):
    """Return True if the stored pipeline PID is still a live process."""
    if not pid:
        return False
    try:
        proc = psutil.Process(int(pid))
        if not proc.is_running() or proc.status() == psutil.STATUS_ZOMBIE:
            return False
        return True
    except (psutil.NoSuchProcess, psutil.AccessDenied, TypeError, ValueError):
        return False


def _clear_run_state():
    process_status["running"] = False
    process_status["pid"] = None
    process_status["exit_code"] = None


def _is_post_pipeline_cmdline(cmd):
    """True for Nextflow/Bakta/work tasks belonging to this post-assembly app."""
    if not cmd:
        return False
    if "post_assembly.py" in cmd or "assembly.py" in cmd:
        return False
    root = str(base_dir)
    if root not in cmd:
        return False
    # Java Nextflow one-jar always embeds main.nf path for this app.
    if "nextflow" in cmd or "/main.nf" in cmd:
        return True
    if f"{root}/work/" in cmd or ".command." in cmd:
        return True
    return any(n in cmd for n in ("bakta_annot.sh", "circlator"))


def _sweep_pipeline_leftovers():
    """Kill leftover Nextflow/work/bakta processes for this app (all process groups)."""
    victims = []
    for proc in psutil.process_iter(["pid", "cmdline"]):
        try:
            if proc.pid == os.getpid():
                continue
            cmd = " ".join(proc.info.get("cmdline") or [])
            if _is_post_pipeline_cmdline(cmd):
                victims.append(proc)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue

    for proc in victims:
        try:
            for child in proc.children(recursive=True):
                try:
                    child.kill()
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
            proc.kill()
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue

    # Second pass — catch stragglers that re-parented during the first pass.
    time.sleep(0.25)
    for proc in psutil.process_iter(["pid", "cmdline"]):
        try:
            if proc.pid == os.getpid():
                continue
            cmd = " ".join(proc.info.get("cmdline") or [])
            if _is_post_pipeline_cmdline(cmd):
                proc.kill()
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue


def _kill_process_tree(pid):
    """Force-kill a shell/Nextflow tree quickly (no long waits)."""
    if pid:
        try:
            pid = int(pid)
        except (TypeError, ValueError):
            pid = None

    if pid:
        # Process-group kill first (started with start_new_session=True).
        for sig in (signal.SIGTERM, signal.SIGKILL):
            try:
                os.killpg(pid, sig)
            except (ProcessLookupError, PermissionError, OSError):
                pass

        try:
            parent = psutil.Process(pid)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            parent = None

        if parent is not None:
            try:
                children = parent.children(recursive=True)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                children = []
            for proc in children + [parent]:
                try:
                    proc.kill()
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass

    _sweep_pipeline_leftovers()
    return True


def run_bact(
    command,
    process_status,
    output_queue,
    output_history,
    work_root=None,
    out_dir=None,
    cleanup_work=False,
):
    """Function to run bactflow"""

    process_status["running"] = True
    history_cap = 1500
    output_history[:] = []
    seen_lines = set()
    stop_event = threading.Event()
    log_thread = None
    if work_root:
        log_thread = threading.Thread(
            target=_tail_nextflow_log,
            args=(work_root, output_queue, stop_event, seen_lines),
            daemon=True,
        )
        log_thread.start()

    master_fd, slave_fd = pty.openpty()
    try:
        process = subprocess.Popen(
            command,
            shell=True,
            executable="/bin/bash",
            stdout=slave_fd,
            stderr=slave_fd,
            stdin=slave_fd,
            close_fds=True,
            start_new_session=True,
        )
    finally:
        os.close(slave_fd)

    process_status["pid"] = process.pid
    buffer = ""

    while True:
        still_marked_running = True
        try:
            still_marked_running = bool(process_status["running"])
        except Exception:
            still_marked_running = True

        if (not still_marked_running) and process.poll() is None:
            _kill_process_tree(process.pid)
            try:
                process.wait(timeout=2)
            except Exception:
                try:
                    process.kill()
                except Exception:
                    pass
            break

        if process.poll() is not None:
            try:
                while True:
                    chunk = os.read(master_fd, 4096).decode(errors="replace")
                    if not chunk:
                        break
                    buffer += chunk
            except OSError:
                pass
            break

        ready, _, _ = select.select([master_fd], [], [], 0.2)
        if master_fd in ready:
            try:
                chunk = os.read(master_fd, 4096).decode(errors="replace")
            except OSError:
                chunk = ""
            if not chunk:
                break
            buffer += chunk

        while buffer:
            split_at = len(buffer)
            for sep in ("\n", "\r"):
                idx = buffer.find(sep)
                if idx != -1:
                    split_at = min(split_at, idx)
            if split_at == len(buffer):
                break
            line = clean_log_line(buffer[:split_at])
            buffer = buffer[split_at + 1 :]
            if not line:
                continue
            output_queue.put(line)
            output_history.append(line)
            if len(output_history) >= history_cap + 250:
                del output_history[0:250]

    if buffer.strip():
        line = clean_log_line(buffer)
        if line:
            output_queue.put(line)
            output_history.append(line)

    os.close(master_fd)
    exit_code = process.wait()
    stop_event.set()
    if log_thread:
        log_thread.join(timeout=2)

    process_status["exit_code"] = int(exit_code) if exit_code is not None else 1
    if exit_code == 0:
        output_queue.put("BactFlow: Nextflow finished successfully.")
        if cleanup_work and out_dir:
            try:
                from work_cleanup import cleanup_run_temps

                output_queue.put(cleanup_run_temps(out_dir))
            except Exception as exc:
                output_queue.put(f"BactFlow: cleanup failed: {exc}")
    else:
        output_queue.put(f"BactFlow: Nextflow exited with code {exit_code}.")

    process_status["running"] = False
    process_status["pid"] = None
    try:
        _lock_release(pid=os.getpid())
    except Exception:
        pass
    return exit_code


# run bactflow

@app.route("/validate-paths", methods=["POST"])
def validate_paths_route():
    """Preflight: genome / DB / output directories must exist (and hold required files)."""
    from path_validate import validate_form_paths

    check = validate_form_paths("post", request.form)
    return jsonify(check), (200 if check.get("ok") else 400)


@app.route('/run_bactflow', methods=['POST', 'GET'])
def run_bactflow():
    """Startp the BactFlow process"""
    if request.method == "POST":
        action = request.args.get("action-assem")
        command = None

        if  action == "run":
            pid = process_status.get("pid")
            if process_status.get("running") and _pipeline_alive(pid):
                return "Bactflow is already running!\n", 400
            if process_status.get("running") or pid:
                _clear_run_state()
            
            
            setup_only = request.form.get("setup_only", 'false')
            cpus = request.form.get('cpus', 1)   
            out_dir = os.path.abspath((request.form.get('out_dir') or './bactflow_out').strip())
            genome_extension = request.form.get('genome_extension', 'fasta')
            checkm_lineag_check = request.form.get('checkm_lineag_check', 'false')
            circle_genome = request.form.get('circle_genome', 'false')  
            tax_class  = request.form.get('tax_class', 'false')
            checkm_db = request.form.get('checkm_db', "")
            gtdbtk_data_path = request.form.get('gtdbtk_data_path', "")
            genome_dir = (request.form.get("genome_dir") or "").strip()
            from path_validate import validate_form_paths

            preflight = validate_form_paths("post", request.form)
            if not preflight.get("ok"):
                return (preflight.get("message") or "Invalid directory path.") + "\n", 400
            if str(setup_only).lower() != "true":
                genome_dir = os.path.abspath(genome_dir)
            run_checkm = request.form.get("run_checkm", "false")
            run_quast = request.form.get("run_quast", "true")
            resume_run = request.form.get("resume_run", "false")
            cleanup_work_dir = request.form.get("cleanup_work_dir", "false")
            run_bakta = request.form.get("bakta_annot", "false")
            bakta_db = request.form.get("bakta_data_path", "")
            run_plasmids = request.form.get("run_plasmids", "false")
            genomad_db = (request.form.get("genomad_db") or "").strip()
            run_mlst = request.form.get("run_mlst", "false")
            command = None
            nf_work = _nf_work_dir(out_dir)

            command = f"""mkdir -p '{out_dir}' '{nf_work}'
                    export NXF_WORK='{nf_work}'
                    cd '{base_dir}' && \\
                    nextflow run {base_dir}/main.nf \\
                    -w '{nf_work}' \\
                    --setup_only {setup_only} \\
                    --cpus {cpus} \\
                    --out_dir '{out_dir}' \\
                    --genome_extension {genome_extension} \\
                    --checkm_lineag_check {checkm_lineag_check} \\
                    --run_flye false \\
                    --circle_genome {circle_genome} \\
                    --tax_class {tax_class} \\
                    --run_checkm {run_checkm} \\
                    --checkm_db '{checkm_db}' \\
                    --bakta_annot {run_bakta} \\
                    --bakta_db '{bakta_db}' \\
                    --gtdbtk_data_path '{gtdbtk_data_path}' \\
                    --run_quast {run_quast} \\
                    --run_plasmids {run_plasmids} \\
                    --genomad_db '{genomad_db}' \\
                    --run_mlst {run_mlst} \\
                    --genome_dir '{genome_dir}' \\
                    -ansi-log true"""
            if _truthy(resume_run, False):
                command = command + " -resume"
            mem_gb, mem_capped, mem_host, _max_gb = _memory_from_form()
            from memory_cap import inject_nextflow_memory, preamble
            command = inject_nextflow_memory(command, out_dir, cpus, mem_gb)
            command = preamble(mem_gb, mem_host, mem_capped) + command
            command = with_nextflow_java(command, mem_gb)
            ok, err = _lock_acquire("Post-assembly")
            if not ok:
                return err + "\n", 409
                
            output_history[:] = []
            process_status["running"] = True
            back_process = Process(
                target=run_bact,
                args=(
                    command,
                    process_status,
                    output_queue,
                    output_history,
                    base_dir,
                    out_dir,
                    _truthy(cleanup_work_dir, False),
                ),
            )
            _start_worker_process(back_process)
            _lock_worker(back_process.pid)
        
            command = None
            return "Bactflow started successfully!\n", 200
            
        if action == "help":
            command = with_nextflow_java(
                f"cd '{base_dir}' && nextflow run {base_dir}/main.nf --help -ansi-log true"
            )
            ok, err = _lock_acquire("Post-assembly")
            if not ok:
                return err + "\n", 409
            output_history[:] = []
            process_status["running"] = True
            back_process = Process(target=run_bact, args=(command, process_status, output_queue, output_history, base_dir))
            _start_worker_process(back_process)
            _lock_worker(back_process.pid)
            
            command = None
            return "Bactflow started successfully!\n", 200
        
        if action == "stop":
            pid = process_status.get("pid")
            was_running = bool(process_status.get("running")) or _pipeline_alive(pid)
            _clear_run_state()
            _lock_release()
            try:
                if pid:
                    _kill_process_tree(pid)
                for proc in psutil.process_iter(["pid", "cmdline"]):
                    try:
                        cmd = " ".join(proc.info.get("cmdline") or [])
                        if str(base_dir) in cmd and (
                            "nextflow" in cmd or "main.nf" in cmd or "bakta_annot.sh" in cmd
                        ):
                            proc.kill()
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        continue
                try:
                    output_queue.put("BactFlow: stopped by user.")
                except Exception:
                    pass
                if was_running or pid:
                    return "Bactflow stopped successfully!\n", 200
                return "No running process to stop.\n", 200
            except Exception as e:
                _clear_run_state()
                return f"Error stopping process: {str(e)}\n", 500

        return f"Unknown action: {action}\n", 400

# a constant ouput
@app.route('/bactflow_output', methods = ['GET'])
def bactflow_output():
    """Returns a constant output"""
    
    output_list = list(output_history)
    return jsonify({"output": output_list})

# Bactflow running status
@app.route('/bactflow_status', methods = ['GET'])
def bactflow_status():
    payload = current_sys_stats()
    payload["status"] = "running" if payload.get("running") else "stopped"
    return jsonify(payload), 200


def _manager_get(key, default=None):
    try:
        return process_status[key]
    except Exception:
        return default


def _job_resource_stats(pid):
    """CPU cores and RSS for the Nextflow process tree."""
    global _JOB_CPU_CACHE
    try:
        parent = psutil.Process(int(pid))
        procs = [parent] + parent.children(recursive=True)
    except (psutil.NoSuchProcess, psutil.AccessDenied, TypeError, ValueError):
        _JOB_CPU_CACHE = {"pid": None, "stamp": 0.0, "cpu_sec": 0.0}
        return None

    cpu_sec = 0.0
    rss = 0
    for proc in procs:
        try:
            times = proc.cpu_times()
            cpu_sec += times.user + times.system
            rss += proc.memory_info().rss
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue

    now = time.time()
    cores = 0.0
    if _JOB_CPU_CACHE["pid"] == pid and _JOB_CPU_CACHE["stamp"]:
        dt = max(now - _JOB_CPU_CACHE["stamp"], 1e-3)
        cores = max(cpu_sec - _JOB_CPU_CACHE["cpu_sec"], 0.0) / dt
    _JOB_CPU_CACHE = {"pid": pid, "stamp": now, "cpu_sec": cpu_sec}
    ncpu = psutil.cpu_count() or 1
    return {
        "job_cores": round(cores, 2),
        "job_cpu_percent": round(min(100.0, cores / ncpu * 100.0), 1),
        "job_rss_gb": round(rss / (1024 ** 3), 2),
        "job_procs": len(procs),
    }


def _cgroup_memory_bases():
    bases = []
    try:
        with open("/proc/self/cgroup", "r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if line.startswith("0::"):
                    rel = line.split("::", 1)[1]
                    bases.append("/sys/fs/cgroup" + (rel if rel.startswith("/") else "/" + rel))
                elif ":memory:" in line:
                    rel = line.split(":")[-1]
                    bases.append("/sys/fs/cgroup/memory" + (rel if rel.startswith("/") else "/" + rel))
    except OSError:
        pass
    bases.extend(["/sys/fs/cgroup", "/sys/fs/cgroup/memory"])
    seen = set()
    out = []
    for base in bases:
        base = base.rstrip("/") or "/"
        if base not in seen:
            seen.add(base)
            out.append(base)
    return out


def _parse_memory_size(text):
    if text is None:
        return None
    raw = str(text).strip().lower()
    if not raw or raw.startswith("auto") or raw in ("max", "unlimited", "inf"):
        return None
    mult = 1
    if raw[-1] in "bkmgt":
        unit = raw[-1]
        num = raw[:-1]
        mult = {"b": 1, "k": 1024, "m": 1024 ** 2, "g": 1024 ** 3, "t": 1024 ** 4}[unit]
    elif raw.endswith(("kb", "mb", "gb", "tb", "ki", "mi", "gi", "ti")):
        unit = raw[-2:]
        num = raw[:-2]
        mult = {
            "kb": 1000, "mb": 1000 ** 2, "gb": 1000 ** 3, "tb": 1000 ** 4,
            "ki": 1024, "mi": 1024 ** 2, "gi": 1024 ** 3, "ti": 1024 ** 4,
        }[unit]
    else:
        num = raw
    try:
        return int(float(num) * mult)
    except ValueError:
        return None


def _env_docker_memory_limit():
    for key in ("BACTFLOW_DOCKER_MEMORY", "BACTFLOW_MEMORY"):
        limit = _parse_memory_size(os.environ.get(key))
        if limit and limit > 0:
            return limit
    return None


def _cgroup_memory_bytes():
    host_total = None
    try:
        with open("/proc/meminfo", "r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("MemTotal:"):
                    host_total = int(line.split()[1]) * 1024
                    break
    except (OSError, ValueError):
        pass

    for base in _cgroup_memory_bases():
        for used_name, lim_name in (
            ("memory.current", "memory.max"),
            ("memory.usage_in_bytes", "memory.limit_in_bytes"),
        ):
            used_path = os.path.join(base, used_name)
            lim_path = os.path.join(base, lim_name)
            try:
                with open(used_path, "r", encoding="utf-8") as handle:
                    used = int(handle.read().strip())
                with open(lim_path, "r", encoding="utf-8") as handle:
                    raw = handle.read().strip()
                if raw in ("max", ""):
                    continue
                limit = int(raw)
                if limit <= 0 or limit >= (1 << 62):
                    continue
                if host_total and limit >= host_total * 0.98:
                    continue
                return used, limit
            except (OSError, ValueError):
                continue
    return None, None


def _docker_aware_memory(fallback_used, fallback_total, fallback_percent):
    used_b, limit_b = _cgroup_memory_bytes()
    env_limit = _env_docker_memory_limit()
    if limit_b is None and env_limit:
        limit_b = env_limit
        if used_b is None:
            used_b = min(int(fallback_used * (1024 ** 3)), limit_b) if fallback_total else 0
    if used_b is not None and limit_b is not None and limit_b > 0:
        return (
            round(used_b / (1024 ** 3), 2),
            round(limit_b / (1024 ** 3), 2),
            round(min(100.0, used_b / limit_b * 100.0), 1),
            "docker",
        )
    return fallback_used, fallback_total, fallback_percent, "host"


def _build_sys_stats():
    vm = psutil.virtual_memory()
    ncpu = psutil.cpu_count(logical=True) or 1
    cpu = round(psutil.cpu_percent(interval=None), 1)
    ram_percent = round(vm.percent, 1)
    ram_used_gb = round(vm.used / (1024 ** 3), 2)
    ram_total_gb = round(vm.total / (1024 ** 3), 2)
    ram_used_gb, ram_total_gb, ram_percent, mem_src = _docker_aware_memory(
        ram_used_gb, ram_total_gb, ram_percent
    )
    payload = {
        "cpu_percent": cpu,
        "cpu_count": ncpu,
        "ram_percent": ram_percent,
        "ram_used_gb": ram_used_gb,
        "ram_total_gb": ram_total_gb,
        "ram_scope": mem_src,
        "running": bool(_manager_get("running", False)),
        "job_cores": None,
        "job_cpu_percent": None,
        "job_rss_gb": None,
        "job_procs": 0,
    }
    pid = _manager_get("pid")
    if payload["running"] and pid:
        job = _job_resource_stats(pid)
        if job:
            payload.update(job)
    return payload


_SYS_STATS = {}
_SYS_STATS_LOCK = threading.Lock()
_SYS_STATS_JSON = os.path.join(base_dir, "static", "sys_stats.json")


def _write_sys_stats_file(payload):
    try:
        os.makedirs(os.path.dirname(_SYS_STATS_JSON), exist_ok=True)
        tmp = _SYS_STATS_JSON + ".tmp"
        with open(tmp, "w", encoding="utf-8") as handle:
            json.dump(payload, handle)
        os.replace(tmp, _SYS_STATS_JSON)
    except OSError:
        pass


def current_sys_stats():
    with _SYS_STATS_LOCK:
        if _SYS_STATS:
            return dict(_SYS_STATS)
    try:
        return _build_sys_stats()
    except Exception as exc:
        return {
            "cpu_percent": 0,
            "cpu_count": 1,
            "ram_percent": 0,
            "ram_used_gb": 0,
            "ram_total_gb": 0,
            "running": False,
            "error": str(exc),
        }


def _sys_stats_loop():
    psutil.cpu_percent(interval=None)
    while True:
        time.sleep(1)
        try:
            payload = _build_sys_stats()
            with _SYS_STATS_LOCK:
                _SYS_STATS.clear()
                _SYS_STATS.update(payload)
            _write_sys_stats_file(payload)
        except Exception:
            continue


threading.Thread(target=_sys_stats_loop, daemon=True, name="sys-stats").start()


@app.route("/sys_stats", methods=["GET"])
def sys_stats():
    return jsonify(current_sys_stats())

#Stream    
@app.route('/stream_bactflow', methods = ['POST', 'GET'])
def stream_bactflow():

  
    def generate():
        last_stats = 0
        # fitst show the history
        for line in output_history:
            yield f"data: {line}\n\n"

        # now stream new output
        while process_status['running']  or not output_queue.empty():
            now = time.time()
            if now - last_stats >= 1:
                yield f"event: stats\ndata: {json.dumps(current_sys_stats())}\n\n"
                last_stats = now
            try:
                line = output_queue.get(timeout=0.5)# wait for output
                yield f"data: {line.strip()}\n\n"
            except Exception:
                if not process_status['running']:
                    break
                continue

        # Do not claim success when Nextflow (or the worker) failed — the UI keeps
        # the terminal open for debugging and must not mark stages as done.
        exit_code = process_status.get("exit_code")
        if exit_code is None:
            exit_code = 0 if not process_status.get("running") else 1
        try:
            exit_code = int(exit_code)
        except (TypeError, ValueError):
            exit_code = 1
        if exit_code == 0:
            yield "data: Process completed\n\n"
        else:
            yield f"data: Process failed (exit {exit_code})\n\n"
 

    return Response(generate(), content_type='text/event-stream')

# progress bar
prog = {"completed": 0}

@app.route("/progress", methods = ["POST", "GET"])
def progress():
    
    if request.method == "POST":
        out_dir = request.form.get('out_dir', './bactflow_out')



        pars = {
            "circle_genome": True if request.form.get("circle_genome") == "true" else False,
            "run_quast": True if request.form.get("run_quast") == "true" else False,
            "tax_class": True if request.form.get("tax_class") == "true" else False,
            "bakta_annot": True if request.form.get("bakta_annot") == "true" else False,
            "run_checkm": True if request.form.get("run_checkm") == "true" else False,
            "run_plasmids": True if request.form.get("run_plasmids") == "true" else False,
            "run_mlst": True if request.form.get("run_mlst") == "true" else False,
            }

       
        def out_dir_watcher(out_dir):
            
            global prog
            
            expected_counts = sum(value is True for value in pars.values())

            pr_stat = process_status["running"]

            if expected_counts > 0:
                while pr_stat:
                    try:
                        current_count = len(os.listdir(out_dir)) if os.path.exists(out_dir) else 0
                        prog = 0
                        prog = {"completed": f"{current_count/expected_counts*100}"}

                        if current_count >= expected_counts:
                            pr_stat = False
                    except Exception as e:
                        print(f"Error watching the output directory: {e}")
                    time.sleep(2)
        wather_thread = threading.Thread(target=out_dir_watcher, args=(out_dir,), daemon=True)
        wather_thread.start() 

    return jsonify(prog) 

# quast — report.html needs sibling JS/CSS; serve the folder under a stable URL.
_QUAST_FS = {}


def _find_quast_dir(out_dir):
    if not out_dir:
        return None
    path = os.path.join(out_dir, "quast_stat")
    if os.path.isdir(path) and os.path.isfile(os.path.join(path, "report.html")):
        return os.path.abspath(path)
    return None


def _quast_mount(out_dir):
    quast_path = _find_quast_dir(out_dir)
    if not quast_path:
        return None, None
    key = hashlib.sha1(quast_path.encode("utf-8")).hexdigest()[:16]
    _QUAST_FS[key] = quast_path
    return key, quast_path


@app.route("/check-quast", methods = ["POST"])  
def check_quast():
    out_dir = request.form.get("out_dir")
    key, quast_path = _quast_mount(out_dir)
    if not quast_path:
        return jsonify({"exists": False})
    return jsonify({
        "exists": True,
        "report_url": url_for("quast_fs", key=key, filename="report.html"),
        "contig_url": url_for("quast_fs", key=key, filename="icarus_viewers/contig_size_viewer.html"),
    })


@app.route("/quast-fs/<key>/<path:filename>")
def quast_fs(key, filename):
    root = _QUAST_FS.get(key)
    if not root or not os.path.isdir(root):
        abort(404)
    root = os.path.abspath(root)
    target = os.path.abspath(os.path.join(root, filename))
    if not (target == root or target.startswith(root + os.sep)):
        abort(403)
    if not os.path.isfile(target):
        abort(404)
    return send_from_directory(root, filename)


@app.route("/check-circ", methods=["POST"])
def check_circ():
    out_dir = (request.form.get("out_dir") or "").strip()
    if not out_dir:
        return jsonify({"exists": False, "files": []})
    circ_dir = os.path.join(out_dir, "circulated_fasta")
    files = []
    if os.path.isdir(circ_dir):
        for name in sorted(os.listdir(circ_dir)):
            low = name.lower()
            if low.endswith((".fasta", ".fa", ".fna")):
                files.append(name)
    return jsonify({
        "exists": bool(files),
        "files": files,
        "dir": circ_dir if files else None,
        "count": len(files),
    })


@app.route("/quast-report", methods = ["POST"])
def quast_report():
    out_dir = request.form.get("out_dir")
    key, quast_path = _quast_mount(out_dir)
    if not quast_path:
        return jsonify({"error": "QUAST report does not exist"}), 404
    return send_from_directory(quast_path, "report.html")

@app.route("/contig-report", methods = ["POST"])
def contig_report():
    out_dir = request.form.get("out_dir")
    key, quast_path = _quast_mount(out_dir)
    if quast_path:
        contig = os.path.join(quast_path, "icarus_viewers/contig_size_viewer.html")
        if os.path.isfile(contig):
            return send_from_directory(quast_path,"icarus_viewers/contig_size_viewer.html")
    return jsonify({"error": "Contig size report does not exist"}), 404

# Bakta
@app.route("/check-bakta", methods = ["POST"])  
def check_bakta():
    out_dir = request.form.get("out_dir")
    gt = request.form.getlist("gene_type")
    gene_type = ",".join(gt)

    if not out_dir:
        return jsonify({"exists": False})
    try:
        gene_count = _ensure_gene_count(out_dir, gene_type)
        if not gene_count:
            return jsonify({"exists": False})

        df = pd.read_csv(gene_count, sep="\t")
        if df.empty:
            return jsonify({"exists": False})
        table_data = df.to_dict(orient="records")
        table_html = """
                <table id="{{ id }}" class="display table table-striped table-bordered nowrap table-hover">
                        <thead>
                            <tr>{% for column in table[0].keys() %}<th>{{ column }}</th>{% endfor %}</tr>
                        </thead>
                        <tbody>
                            {% for row in table %}
                            <tr>{% for value in row.values() %}<td>{{ value }}</td>{% endfor %}</tr>
                            {% endfor %}
                        </tbody>
                    </table>
                """
        return jsonify({
            "exists": True,
            "count_tab": render_template_string(
                table_html,
                id="bakta-tab",
                tabnumber="Table 1: Count table of annotated genes by Bakta.",
                tabcaption="This is the abundance of all genes (selected based on gene type input) at all genomes.",
                table=table_data,
            ),
        })
    except subprocess.CalledProcessError as e:
        return jsonify({
            "exists": False,
            "error": f"Gene counter failed with error code {e.returncode}",
            "stderr": (e.stderr or "").strip(),
            "stdout": (e.stdout or "").strip(),
        })
    except Exception:
        return jsonify({"exists": False})

def _find_suffix_files(root, suffixes):
    found = []
    if not root or not os.path.isdir(root):
        return found
    for dirpath, _dirnames, names in os.walk(root):
        for name in names:
            low = name.lower()
            if any(low.endswith(sfx) for sfx in suffixes):
                found.append(os.path.join(dirpath, name))
    return found


def _ensure_gene_count(out_dir, gene_type):
    """Build gene_count.tsv from Bakta TSVs when missing. Returns path or None."""
    if not out_dir:
        return None
    bakta_path = os.path.join(out_dir, "bakta_out")
    gene_count = os.path.join(out_dir, "gene_count.tsv")
    if not os.path.isdir(bakta_path):
        return None
    if os.path.isfile(gene_count) and os.path.getsize(gene_count) > 0:
        return gene_count
    tsvs = _find_suffix_files(bakta_path, (".tsv",))
    if not tsvs:
        return None
    types = gene_type or "cds"
    command = f"""
            cp {bakta_path}/*/*.tsv {bakta_path}
            rm {bakta_path}/*inference.tsv {bakta_path}/*hypotheticals.tsv
            mkdir -p {bakta_path}/genes && mv {bakta_path}/*.tsv {bakta_path}/genes

            python {base_dir}/gene_counter_bakta.py -d {bakta_path}/genes -t {types} -o {out_dir}/gene_count
            """
    subprocess.run(command, shell=True, text=True, check=True, capture_output=True)
    if os.path.isfile(gene_count) and os.path.getsize(gene_count) > 0:
        return gene_count
    return None


_HEATMAP_HYPO_RE = re.compile(
    r"hypothetical|uncharacterized|unknown function|putative protein|unannotated",
    re.I,
)


def _gene_heatmap_payload(gene_count_path, n=50):
    df = pd.read_csv(gene_count_path, sep="\t")
    if df.empty or df.shape[1] < 2:
        return None
    product_col = df.columns[0]
    genomes = [col for col in df.columns[1:]]
    numeric = df[genomes].apply(pd.to_numeric, errors="coerce").fillna(0)
    products = df[product_col].fillna("").astype(str).str.strip()
    named = ~products.str.contains(_HEATMAP_HYPO_RE, na=False)
    keep = named & (products != "") & (numeric.sum(axis=1) > 0)
    numeric = numeric.loc[keep]
    products = products.loc[keep]
    if numeric.empty:
        return None
    rank = numeric.sum(axis=1) if len(genomes) < 2 else numeric.var(axis=1)
    raw_n = str(n).strip().lower() if n is not None else "50"
    total = int(len(numeric))
    if raw_n in {"all", "0", "-1"}:
        take = total
    else:
        try:
            take = int(float(raw_n))
        except (TypeError, ValueError):
            take = 50
        take = max(1, min(take, total))
    n = take
    top_idx = rank.sort_values(ascending=False).head(n).index
    sub = numeric.loc[top_idx]
    return {
        "exists": True,
        "n": n,
        "n_total": int(len(numeric)),
        "n_genomes": len(genomes),
        "rank": "abundance" if len(genomes) < 2 else "variance",
        "genomes": genomes,
        "genes": products.loc[top_idx].tolist(),
        "abundance": sub.round(0).astype(int).values.tolist(),
        "prevalence": (sub > 0).astype(int).values.tolist(),
        "variance": [round(float(v), 4) for v in rank.loc[top_idx].tolist()],
    }


@app.route("/heatmap-genes", methods=["POST"])
def heatmap_genes():
    out_dir = request.form.get("out_dir", "").strip()
    gt = request.form.getlist("gene_type")
    gene_type = ",".join(gt)
    n = request.form.get("n") or request.form.get("heatmap_n") or 50
    if not out_dir:
        return jsonify({"exists": False, "error": "Set an output directory first."})
    try:
        gene_count = _ensure_gene_count(out_dir, gene_type)
        if not gene_count:
            return jsonify({"exists": False})
        payload = _gene_heatmap_payload(gene_count, n)
        if not payload:
            return jsonify({"exists": False, "error": "Gene count table is empty."})
        return jsonify(payload)
    except subprocess.CalledProcessError as e:
        return jsonify({
            "exists": False,
            "error": f"Gene counter failed with error code {e.returncode}",
            "stderr": (e.stderr or "").strip(),
        })
    except Exception as exc:
        return jsonify({"exists": False, "error": str(exc) or "Could not build heatmap."})


@app.route("/check-bakta-ready", methods=["POST"])
def check_bakta_ready():
    """Check Bakta outputs used by the circular plot (gbk/gbff/gff)."""
    out_dir = request.form.get("out_dir", "").strip()
    if not out_dir:
        return jsonify({
            "ready": False,
            "plot_ready": False,
            "message": "Set an output directory first.",
            "plot_message": "Set an output directory first.",
        })

    bakta_path = os.path.join(out_dir, "bakta_out")
    if not os.path.isdir(bakta_path):
        msg = "No Bakta output yet. Run BactFlow with Bakta enabled."
        return jsonify({
            "ready": False,
            "plot_ready": False,
            "message": msg,
            "plot_message": msg,
        })

    plot_files = _find_suffix_files(bakta_path, (".gbff", ".gbk", ".gb", ".gff3", ".gff"))
    plot_ready = len(plot_files) > 0
    kinds = sorted({os.path.splitext(p)[1].lower() for p in plot_files})
    return jsonify({
        "ready": plot_ready,
        "plot_ready": plot_ready,
        "plot_count": len(plot_files),
        "gbk_count": len(plot_files),
        "plot_kinds": kinds,
        "plot_message": None if plot_ready else (
            "Bakta output exists but no .gbk/.gbff/.gff/.gff3 files were found."
        ),
        "message": None if plot_ready else (
            "Bakta output exists but no .gbk/.gbff/.gff/.gff3 files were found."
        ),
    })

# circular
@app.route("/circular", methods=["POST"])
def circular():
    out_dir = request.form.get("out_dir")
    if not out_dir:
        return jsonify({"plot": False, "reason": "missing_out_dir"}), 200
    generate = str(request.form.get("generate", "false")).lower() == "true"
    gbk_dir = os.path.join(out_dir, "bakta_out")
    crc_plt = os.path.join(out_dir, "circular_plot.png")
    params_file = os.path.join(out_dir, "circular_plot_params.json")  
    
    
    params = {
        "add_gc": request.form.get("add_gc"),
        "add_skew": request.form.get("add_skew"),
        "dpi": int(request.form.get("dpi", 200)),
        "figsize": int(request.form.get("figsize", 10)),
        "interval": int(request.form.get("interval", 3)),
        "f_color": request.form.get("f_color", "#1E90FF"),
        "r_color": request.form.get("r_color", "#FF7261"),
        "feature_types": ",".join(_selected_gene_types()),
    }

    last_params = {}


    if os.path.exists(params_file):
        try:
            with open(params_file, 'r') as f:
                if os.path.getsize(params_file) > 0:
                    last_params = json.load(f)

        except (json.JSONDecodeError, IOError) as e:
            print(f"Error reading params file for the plot: {e}")
    
    needs_regen = (not os.path.exists(crc_plt)) or params != last_params or generate
    if generate and os.path.exists(crc_plt):
        try:
            os.remove(crc_plt)
        except OSError:
            pass
        json_old = os.path.join(out_dir, "circular_plot.json")
        if os.path.exists(json_old):
            try:
                os.remove(json_old)
            except OSError:
                pass
        needs_regen = True

    if needs_regen and not generate:
        return jsonify({"plot": False, "reason": "not_generated"}), 200

    if needs_regen:
        print("🔄 Parameters changed or plot missing. Regenerating plot...")
        if not os.path.isdir(gbk_dir):
            return jsonify({"plot": False, "reason": "no_bakta_dir"}), 200
        try:
            with open(params_file, "w") as f:
                json.dump(params, f, indent=4)
        except Exception as e:
            return jsonify({"plot": False, "error": f"Failed to write params: {e}"}), 500

       
        command = f"""
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bactflow
python3 {shlex.quote(os.path.join(base_dir, "circular_plotter.py"))} -d {shlex.quote(gbk_dir)} -o {shlex.quote(out_dir)} \
{"--add_gc" if params["add_gc"] == "True" else ""} \
{"--add_skew" if params["add_skew"] == "True" else ""} \
--dpi {params["dpi"]} \
--figsize {params["figsize"]} \
--interval {params["interval"]} \
--f_color {shlex.quote(params["f_color"])} \
--r_color {shlex.quote(params["r_color"])} \
--feature_types {shlex.quote(params.get("feature_types") or "cds")}
"""

        try:
            completed = _run_bash(command)
        except Exception as e:
            return jsonify({"plot": False, "error": str(e)}), 200
        if completed.returncode != 0:
            err = (completed.stderr or completed.stdout or "Circular plotter failed.").strip()
            return jsonify({"plot": False, "error": err[-3000:]}), 200

    if os.path.exists(crc_plt):
        with open(crc_plt, "rb") as img_file:
            img_base = base64.b64encode(img_file.read()).decode("utf-8")
        payload = {"plot": f"data:image/png;base64,{img_base}"}
        return jsonify(payload)
    else:
        return jsonify({"plot": False, "reason": "not_found_after_generation"}), 200


@app.route("/linear", methods=["POST"])
def linear():
    out_dir = request.form.get("out_dir")
    if not out_dir:
        return jsonify({"plot": False, "reason": "missing_out_dir"}), 200
    gbk_dir = os.path.join(out_dir, "bakta_out")
    if not os.path.isdir(gbk_dir):
        return jsonify({"plot": False, "reason": "no_bakta_dir"}), 200
    try:
        from linear_plotter import build_linear_figure
        figure, meta = build_linear_figure(
            gbk_dir,
            feature_types=_selected_gene_types(),
            genome=(request.form.get("linear_genome") or "").strip(),
            contig=(request.form.get("linear_contig") or "").strip(),
        )
    except Exception as exc:
        return jsonify({"plot": False, "error": str(exc)}), 200
    if not figure:
        return jsonify({"plot": False, **(meta or {})}), 200
    payload = {"plot": True, **(meta or {})}
    if not (isinstance(figure, dict) and figure.get("canvas")):
        payload["figure"] = figure
    return jsonify(payload), 200


@app.route("/linear-links", methods=["POST"])
def linear_links():
    out_dir = request.form.get("out_dir")
    if not out_dir:
        return jsonify({"links": [], "error": "missing_out_dir"})
    gbk_dir = os.path.join(out_dir, "bakta_out")
    if not os.path.isdir(gbk_dir):
        return jsonify({"links": [], "error": "no_bakta_dir"})
    try:
        from linear_plotter import build_synteny_links
        payload = build_synteny_links(
            gbk_dir,
            genome=(request.form.get("linear_genome") or "").strip(),
            contig=(request.form.get("linear_contig") or "").strip(),
            mode=(request.form.get("cord_mode") or "pairwise").strip(),
            reference=(request.form.get("cord_ref") or "").strip(),
            min_pid=request.form.get("min_pid") or 30,
            feature_types=_selected_gene_types(),
        )
        return jsonify(payload)
    except Exception as exc:
        return jsonify({"links": [], "error": str(exc)})


@app.route("/linear-seq", methods=["POST"])
def linear_seq():
    out_dir = request.form.get("out_dir")
    if not out_dir:
        return jsonify({"seq": ""})
    gbk_dir = os.path.join(out_dir, "bakta_out")
    if not os.path.isdir(gbk_dir):
        return jsonify({"seq": ""})
    try:
        from linear_plotter import slice_sequence
        payload = slice_sequence(
            gbk_dir,
            genome=(request.form.get("linear_genome") or "").strip(),
            contig=(request.form.get("linear_contig") or "").strip(),
            start=request.form.get("start") or 1,
            end=request.form.get("end") or 1,
        )
        return jsonify(payload)
    except Exception as exc:
        return jsonify({"seq": "", "error": str(exc)})


    
@app.route("/taxa-report", methods = ["POST"])
def taxa_report():
    raw_out = (request.form.get("out_dir") or "").strip()
    if not raw_out:
        return jsonify({"exists": False})
    out_dir = os.path.abspath(raw_out)
    try:
        from gtdb_report import abundance_display_rows, load_gtdb_taxonomy, prepare_gtdb_tree, taxonomy_display_rows
        extra_roots = [
            os.path.join(base_dir, "work"),
            os.path.join(base_dir, "bactflow_out", ".nextflow-work"),
            _nf_work_dir(out_dir),
            os.path.join(os.getcwd(), "work"),
            os.path.join(os.environ.get("HOME") or "", "work"),
            os.environ.get("NXF_WORK") or "",
        ]
        df = load_gtdb_taxonomy(out_dir, extra_roots=extra_roots)
        tax_rows = taxonomy_display_rows(df) if df is not None and not df.empty else []
        if not tax_rows:
            return jsonify({"exists": False})
        abund_rows = abundance_display_rows(df)
        table_html = """
            <table id="{{ id }}" class="display table table-striped table-bordered nowrap table-hover">
                <thead>
                    <tr>{% for column in table[0].keys() %}<th>{{ column }}</th>{% endfor %}</tr>
                </thead>
                <tbody>
                    {% for row in table %}
                    <tr>{% for value in row.values() %}<td>{{ value }}</td>{% endfor %}</tr>
                    {% endfor %}
                </tbody>
            </table>
        """
        payload = {
            "exists": True,
            "n_genomes": len(tax_rows),
            "taxa_table": render_template_string(table_html, id="taxa-tab", table=tax_rows),
        }
        if abund_rows:
            payload["abund_table"] = render_template_string(table_html, id="taxa-abund-tab", table=abund_rows)
        newick, tree_source = prepare_gtdb_tree(out_dir, extra_roots=extra_roots)
        payload["has_tree"] = bool(newick)
        payload["newick"] = newick or ""
        payload["tree_source"] = tree_source or ""
        return jsonify(payload)
    except Exception as exc:
        print(f"taxa-report failed: {exc}")
        return jsonify({"exists": False})


@app.route("/check-checkm", methods=["POST"])
def check_checkm():
    out_dir = (request.form.get("out_dir") or "").strip()
    if not out_dir:
        return jsonify({"exists": False, "has_tree": False})
    try:
        from checkm_report import (
            find_checkm_files,
            load_gtdb_species_map,
            parse_checkm_lineage_table,
            prepare_checkm_tree,
            read_text_if_nonempty,
        )
        lineage_path, taxon_path, genome_path = find_checkm_files(out_dir)
        lineage_text = read_text_if_nonempty(lineage_path)
        rows = parse_checkm_lineage_table(lineage_text) if lineage_text else []
        if not rows:
            return jsonify({"exists": False, "has_tree": False})
        table_html = """
            <table id="{{ id }}" class="display table table-striped table-bordered nowrap table-hover">
                <thead>
                    <tr>{% for column in table[0].keys() %}<th>{{ column }}</th>{% endfor %}</tr>
                </thead>
                <tbody>
                    {% for row in table %}
                    <tr>{% for value in row.values() %}<td>{{ value }}</td>{% endfor %}</tr>
                    {% endfor %}
                </tbody>
            </table>
        """
        gtdb_map = load_gtdb_species_map(out_dir)
        newick = ""
        tree_source = ""
        for path, kind in ((taxon_path, "taxon_tree.newick"), (genome_path, "genome_tree")):
            raw = read_text_if_nonempty(path)
            if not raw:
                continue
            prepared = prepare_checkm_tree(raw, rows, gtdb_map)
            if prepared:
                newick = prepared
                tree_source = kind
                break
        return jsonify({
            "exists": True,
            "checkm_table": render_template_string(table_html, id="checkm-tab", table=rows),
            "has_tree": bool(newick),
            "newick": newick,
            "tree_source": tree_source,
            "n_genomes": len(rows),
        })
    except Exception as exc:
        print(f"check-checkm failed: {exc}")
        return jsonify({"exists": False, "has_tree": False})


@app.route("/check-mlst", methods=["POST"])
def check_mlst():
    out_dir = (request.form.get("out_dir") or "").strip()
    if not out_dir:
        return jsonify({"exists": False})
    try:
        from mlst_report import load_mlst_rows, mlst_chart_payload, mlst_display_rows, mlst_run_complete
        abs_out = os.path.abspath(out_dir)
        if not mlst_run_complete(abs_out):
            return jsonify({"exists": False})
        rows = mlst_display_rows(abs_out)
        typed = load_mlst_rows(abs_out)
        table_html = """
            <table id="{{ id }}" class="display table table-striped table-bordered nowrap table-hover">
                <thead>
                    <tr>{% for column in table[0].keys() %}<th>{{ column }}</th>{% endfor %}</tr>
                </thead>
                <tbody>
                    {% for row in table %}
                    <tr>{% for value in row.values() %}<td>{{ value }}</td>{% endfor %}</tr>
                    {% endfor %}
                </tbody>
            </table>
        """
        payload = {
            "exists": True,
            "n_genomes": len(rows),
            "n_typed": sum(1 for r in typed if str(r.get("ST") or "-") not in ("-", "", "0")),
            "mlst_table": "",
            "charts": mlst_chart_payload(typed),
        }
        if rows:
            payload["mlst_table"] = render_template_string(table_html, id="mlst-tab", table=rows)
        return jsonify(payload)
    except Exception as exc:
        print(f"check-mlst failed: {exc}")
        return jsonify({"exists": False})


@app.route("/generate-report", methods=["POST"])
def generate_report():
    out_dir = (request.form.get("out_dir") or "").strip()
    if not out_dir:
        return jsonify({"ok": False, "error": "out_dir is required"}), 400
    abs_out = os.path.abspath(out_dir)
    if not os.path.isdir(abs_out):
        return jsonify({"ok": False, "error": f"Output directory does not exist: {abs_out}"}), 400
    try:
        from genome_report import build_reports
        result = build_reports(abs_out)
        return jsonify(result)
    except Exception as exc:
        print(f"generate-report failed: {exc}")
        return jsonify({"ok": False, "error": str(exc)}), 500


@app.route("/download-report/<kind>", methods=["GET"])
def download_report(kind):
    out_dir = (request.args.get("out_dir") or "").strip()
    if not out_dir:
        return jsonify({"error": "out_dir is required"}), 400
    abs_out = os.path.abspath(out_dir)
    report_dir = os.path.join(abs_out, "bactflow_report")
    if kind == "pdf":
        name = "bactflow_genome_report.pdf"
        mimetype = "application/pdf"
    elif kind == "html":
        name = "bactflow_genome_report.html"
        mimetype = "text/html; charset=utf-8"
    else:
        return jsonify({"error": "kind must be pdf or html"}), 400
    path = os.path.join(report_dir, name)
    if not os.path.isfile(path):
        return jsonify({"error": "Report not found. Click Generate report first."}), 404
    return send_from_directory(
        report_dir,
        name,
        as_attachment=True,
        download_name=name,
        mimetype=mimetype,
    )


@app.route("/check-plasmids", methods=["POST"])
def check_plasmids():
    out_dir = (request.form.get("out_dir") or "").strip()
    if not out_dir:
        return jsonify({"exists": False})
    try:
        from plasmid_report import (
            find_plasmid_summary,
            find_virus_summary,
            genomad_run_complete,
            genomad_scanned_genomes,
            plasmid_chart_payload,
            plasmid_display_rows,
            virus_chart_payload,
            virus_display_rows,
        )
        abs_out = os.path.abspath(out_dir)
        if not genomad_run_complete(abs_out):
            return jsonify({"exists": False})
        rows = plasmid_display_rows(abs_out)
        virus_rows = virus_display_rows(abs_out)
        genomes = genomad_scanned_genomes(abs_out)
        has_plasmid_summary = bool(find_plasmid_summary(abs_out))
        has_virus_summary = bool(find_virus_summary(abs_out))
        table_html = """
            <table id="{{ id }}" class="display table table-striped table-bordered nowrap table-hover">
                <thead>
                    <tr>{% for column in table[0].keys() %}<th>{{ column }}</th>{% endfor %}</tr>
                </thead>
                <tbody>
                    {% for row in table %}
                    <tr>{% for value in row.values() %}<td>{{ value }}</td>{% endfor %}</tr>
                    {% endfor %}
                </tbody>
            </table>
        """
        payload = {
            "exists": True,
            "n_plasmids": len(rows),
            "n_viruses": len(virus_rows),
            "genomes": genomes,
            "plasmid_empty": len(rows) == 0 and (has_plasmid_summary or has_virus_summary),
            "virus_empty": len(virus_rows) == 0 and (has_plasmid_summary or has_virus_summary),
            "plasmid_table": "",
            "virus_table": "",
            "charts": plasmid_chart_payload(rows) if rows else None,
            "virus_charts": virus_chart_payload(virus_rows) if virus_rows else None,
        }
        if rows:
            payload["plasmid_table"] = render_template_string(
                table_html, id="plasmid-tab", table=rows
            )
        if virus_rows:
            payload["virus_table"] = render_template_string(
                table_html, id="virus-tab", table=virus_rows
            )
        return jsonify(payload)
    except Exception as exc:
        print(f"check-plasmids failed: {exc}")
        return jsonify({"exists": False})


@app.route("/snp-finder", methods = ["POST"])
def snp_finder():
    out_dir = request.form.get("out_dir")
    genomes = request.form.get("genome_dir")
    reference = request.form.get("ref_genome")
    cpus = request.form.get("cpus")
    if not reference:
        return jsonify({"exists": False, "error": "Select a reference (wild type) genome first."}), 400
    if not genomes:
        return jsonify({"exists": False, "error": "Genome directory is missing."}), 400
    if not out_dir:
        return jsonify({"exists": False, "error": "Output directory is missing."}), 400

    command = f"""
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bactflow
{shlex.quote(os.path.join(base_dir, "variant_finder.sh"))} -r {shlex.quote(reference)} -g {shlex.quote(genomes)} -o {shlex.quote(os.path.join(out_dir, "snps"))} -t {shlex.quote(str(cpus or 1))}
"""
    try:
        completed = _run_bash(command)
    except Exception as exc:
        return jsonify({"exists": False, "error": str(exc)}), 500

    if completed.returncode != 0:
        err = (completed.stderr or completed.stdout or "SNP finder command failed.").strip()
        return jsonify({"exists": False, "error": err[-2000:]}), 500

    collector = collect_variant_tables or _load_collect_variant_tables()
    if collector is None:
        snp_dir = os.path.join(out_dir, "snps")
        snp_files = glob.glob(os.path.join(snp_dir, "**", "*.vcf"), recursive=True)
        if snp_files:
            return jsonify({
                "exists": True,
                "roots": [snp_dir],
                "genomes": [],
                "message": f"VCF files are in {snp_dir} ({len(snp_files)} file(s)). Restart post-assembly to enable table view.",
            })
        return jsonify({"exists": False, "error": "SNP finder finished but no VCF files were found."})

    payload = collector(out_dir, kind="snps")
    if not payload.get("exists"):
        return jsonify({"exists": False, "error": payload.get("error") or "No SNP VCF files found."})
    return jsonify(payload)


@app.route("/snp-report", methods=["POST"])
def snp_report():
    """Load existing SNP VCFs from out_dir without re-running the finder."""
    out_dir = request.form.get("out_dir")
    if not out_dir:
        return jsonify({"exists": False, "error": "Output directory is missing."}), 400
    collector = collect_variant_tables or _load_collect_variant_tables()
    if collector is None:
        snp_dir = os.path.join(out_dir, "snps")
        return jsonify({
            "exists": os.path.isdir(snp_dir),
            "roots": [snp_dir] if os.path.isdir(snp_dir) else [],
            "genomes": [],
            "error": "vcf_viewer module is not available in this server process. Restart post-assembly.",
        }), 500
    payload = collector(out_dir, kind="snps")
    if not payload.get("exists"):
        return jsonify({"exists": False, "error": payload.get("error") or "No SNP VCF files found."})
    return jsonify(payload)


@app.route("/svs-finder", methods = ["POST"])
def svs_finder():
    out_dir = request.form.get("out_dir")
    genomes = request.form.get("genome_dir")
    reference = request.form.get("ref_genome")
    cpus = request.form.get("cpus")
    if not reference:
        return jsonify({"exists": False, "error": "Select a reference (wild type) genome first."}), 400
    if not genomes:
        return jsonify({"exists": False, "error": "Genome directory is missing."}), 400
    if not out_dir:
        return jsonify({"exists": False, "error": "Output directory is missing."}), 400

    command = f"""
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bactflow
{shlex.quote(os.path.join(base_dir, "vc_medaka.sh"))} -r {shlex.quote(reference)} -g {shlex.quote(genomes)} -o {shlex.quote(os.path.join(out_dir, "vcs"))} -c {shlex.quote(str(cpus or 1))}
"""
    try:
        completed = _run_bash(command)
    except Exception as exc:
        return jsonify({"exists": False, "error": str(exc)}), 500

    if completed.returncode != 0:
        err = (completed.stderr or completed.stdout or "Variant calling command failed.").strip()
        return jsonify({"exists": False, "error": err[-2000:]}), 500

    collector = collect_variant_tables or _load_collect_variant_tables()
    if collector is None:
        vcf_dir = os.path.join(out_dir, "vcs")
        vcf_files = glob.glob(os.path.join(vcf_dir, "**", "*.vcf"), recursive=True)
        if vcf_files:
            return jsonify({
                "exists": True,
                "roots": [vcf_dir],
                "genomes": [],
                "message": f"VCF files are in {vcf_dir} ({len(vcf_files)} file(s)). Restart post-assembly to enable table view.",
            })
        return jsonify({"exists": False, "error": "Variant calling finished but no VCF files were found."})

    payload = collector(out_dir, kind="vcs")
    if not payload.get("exists"):
        return jsonify({"exists": False, "error": payload.get("error") or "No variant-call VCF files found."})
    return jsonify(payload)


@app.route("/svs-report", methods=["POST"])
def svs_report():
    """Load existing variant-call VCFs from out_dir without re-running."""
    out_dir = request.form.get("out_dir")
    if not out_dir:
        return jsonify({"exists": False, "error": "Output directory is missing."}), 400
    collector = collect_variant_tables or _load_collect_variant_tables()
    if collector is None:
        return jsonify({"exists": False, "error": "vcf_viewer module is not available. Restart post-assembly."}), 500
    payload = collector(out_dir, kind="vcs")
    if not payload.get("exists"):
        return jsonify({"exists": False, "error": payload.get("error") or "No variant-call VCF files found."})
    return jsonify(payload)

# abundance
@app.route("/abund-run", methods = ["POST"])
def abund_finder():
    out_dir = request.form.get("out_dir")
    gene_files = request.form.get("gene_files")
    prevalance = "false"
    # Enzyme matching needs CDS products/genes; ignore UI gene_type for strain finder.
    gene_type = "cds"
    width = request.form.get("plot_width") or "10"
    height = request.form.get("plot_height") or "10"

    
    enzyme_file = request.form.get("enzyme_loc")
    if not out_dir:
        return jsonify({"exists": False, "error": "Output directory is missing."}), 400
    if not gene_files:
        return jsonify({"exists": False, "error": "Gene annotation directory is missing."}), 400
    if not enzyme_file:
        return jsonify({"exists": False, "error": "Enzyme file path is missing."}), 400
    try:
        width = str(max(5.0, min(40.0, float(width))))
        height = str(max(5.0, min(40.0, float(height))))
    except (TypeError, ValueError):
        width, height = "12", "12"
    output = os.path.join(out_dir, "strain_finder")
    count_tab = os.path.join(output, "abundance.tsv")
    plot = _first_existing(
        os.path.join(output, "requested_genes_abundance.jpeg"),
        os.path.join(output, "requested_genes_abundance.jpg"),
        os.path.join(output, "requested_genes_abundance.png"),
    )
    print(f"this is widht {width} and heigth {height}")
    command = (
        f"{shlex.quote(os.path.join(base_dir, 'strain_finder.sh'))} "
        f"-d {shlex.quote(gene_files or '')} -g {shlex.quote(gene_type)} "
        f"-f abundance -o {shlex.quote(output)} -e {shlex.quote(enzyme_file or '')} "
        f"-c false -p {prevalance} -w {shlex.quote(str(width))} -l {shlex.quote(str(height))}"
    )

    try:
        completed = subprocess.run(command, shell=True, text=True, check=False, capture_output=True)
    except Exception as exc:
        return jsonify({"exists": False, "error": str(exc)}), 500
    if completed.returncode != 0:
        err = (completed.stderr or completed.stdout or "Abundance strain finder failed.").strip()
        return jsonify({"exists": False, "error": err[-3000:]}), 500

    # Re-resolve plot after the run (jpeg and/or png may appear).
    plot = _first_existing(
        os.path.join(output, "requested_genes_abundance.jpeg"),
        os.path.join(output, "requested_genes_abundance.jpg"),
        os.path.join(output, "requested_genes_abundance.png"),
    )
    if os.path.exists(count_tab) and plot:
        df  = pd.read_csv(count_tab, sep = "\t")
        data = df.to_dict(orient = "records")
        table_html = """
                <table id="{{ id }}" class="display table table-striped table-bordered nowrap table-hover">
                        <thead> 
                            <tr>{% for column in table[0].keys() %}<th>{{ column }}</th>{% endfor %}</tr>
                        </thead>
                        <tbody>
                            {% for row in table %}
                            <tr>{% for value in row.values() %}<td>{{ value }}</td>{% endfor %}</tr>
                            {% endfor %}
                        </tbody>
                    </table>
                """
      
        return jsonify({
                "exists": True, 
                 "abund_table": render_template_string(
                    table_html, id="abund-tab", tabnumber = "Table 3: Abundance table of requested genes.", table=data
                ),
                "plot_abund" : _image_data_uri(plot)
                })
    missing = []
    if not os.path.exists(count_tab):
        missing.append("abundance.tsv")
    if not plot:
        missing.append("requested_genes_abundance.(jpeg|png)")
    return jsonify({"exists": False, "error": f"Missing expected outputs: {', '.join(missing)}"}), 400



# prevalance
@app.route("/prev-run", methods = ["POST"])
def prev_finder():
    out_dir = request.form.get("out_dir")
    gene_files = request.form.get("gene_files")
    prevalance = "true"
    # Enzyme matching needs CDS products/genes; ignore UI gene_type for strain finder.
    gene_type = "cds"
    
    enzyme_file = request.form.get("enzyme_loc")
    if not out_dir:
        return jsonify({"exists": False, "error": "Output directory is missing."}), 400
    if not gene_files:
        return jsonify({"exists": False, "error": "Gene annotation directory is missing."}), 400
    if not enzyme_file:
        return jsonify({"exists": False, "error": "Enzyme file path is missing."}), 400
    output = os.path.join(out_dir, "strain_finder")
    count_tab = os.path.join(output, "prevalance.tsv")
    width = request.form.get("plot_width") or "12"
    height = request.form.get("plot_height") or "12"
    try:
        width = str(max(5.0, min(40.0, float(width))))
        height = str(max(5.0, min(40.0, float(height))))
    except (TypeError, ValueError):
        width, height = "12", "12"

    command = (
        f"{shlex.quote(os.path.join(base_dir, 'strain_finder.sh'))} "
        f"-d {shlex.quote(gene_files or '')} -g {shlex.quote(gene_type)} "
        f"-f prevalance -o {shlex.quote(output)} -e {shlex.quote(enzyme_file or '')} "
        f"-c false -p {prevalance} -w {shlex.quote(str(width))} -l {shlex.quote(str(height))}"
    )

    completed = subprocess.run(command, shell=True, text=True, check=False, capture_output=True)
    if completed.returncode != 0:
        err = (completed.stderr or completed.stdout or "Prevalance strain finder failed.").strip()
        return jsonify({"exists": False, "error": err[-3000:]}), 500

    plot = _first_existing(
        os.path.join(output, "requested_genes_prevalence.jpeg"),
        os.path.join(output, "requested_genes_prevalence.jpg"),
        os.path.join(output, "requested_genes_prevalence.png"),
    )
    if os.path.exists(count_tab) and plot:
        df  = pd.read_csv(count_tab, sep = "\t")
        data = df.to_dict(orient = "records")
        table_html = """
                <table id="{{ id }}" class="display table table-striped table-bordered nowrap table-hover">
                        <thead> 
                            <tr>{% for column in table[0].keys() %}<th>{{ column }}</th>{% endfor %}</tr>
                        </thead>
                        <tbody>
                            {% for row in table %}
                            <tr>{% for value in row.values() %}<td>{{ value }}</td>{% endfor %}</tr>
                            {% endfor %}
                        </tbody>
                    </table>
                """
      
        return jsonify({
                "exists": True, 
                 "prev_table": render_template_string(
                    table_html, id="prev-tab", tabnumber = "Table 4: Prevalance table of requested genes.", table=data
                ),
                "plot_prev" : _image_data_uri(plot)
                })
    missing = []
    if not os.path.exists(count_tab):
        missing.append("prevalance.tsv")
    if not plot:
        missing.append("requested_genes_prevalence.(jpeg|png)")
    return jsonify({"exists": False, "error": f"Missing expected outputs: {', '.join(missing)}"}), 400




if __name__ == '__main__':
    # Docker: open via host browser hook. Local: webbrowser/xdg-open.
    # bactflow.sh sets BACTFLOW_NO_BROWSER=1 and opens the tab itself.
    want_browser = (
        os.environ.get("BACTFLOW_NO_BROWSER") != "1"
        or os.environ.get("BACTFLOW_FORCE_BROWSER") == "1"
    )
    if want_browser:
        Timer(1, open_browser).start()
    # Keep Werkzeug debugger off by default — debug + fork breaks Start BactFlow
    # (BrokenPipeError on stdout flush) and returns HTML error pages to the UI.
    flask_debug = os.environ.get("BACTFLOW_FLASK_DEBUG") == "1"
    app.run(
        debug=flask_debug,
        port=_ui_port(5001),
        host="0.0.0.0",
        use_reloader=False,
        threaded=True,
    )
