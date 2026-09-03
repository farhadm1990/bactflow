#!/usr/bin/env python3
import sys, subprocess, importlib, os, shlex
import pty
import select
from flask import Flask, render_template, request, redirect, Response, send_from_directory, stream_with_context, jsonify, render_template_string
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
app = Flask(__name__, 
            template_folder = os.path.join(base_dir, "templates"),
            static_folder = os.path.join(base_dir, "static"))

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


def with_nextflow_java(command):
    return NF_JAVA_SETUP + "\n" + command

ANSI_RE = re.compile(r"\x1b\[[0-9;?]*[A-Za-z]|\x1b\][^\x07]*\x07")


def _run_bash(script):
    return subprocess.run(
        ["bash", "-lc", script],
        check=False,
        text=True,
        capture_output=True,
    )


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


def open_browser():
    webbrowser.open(url="http://127.17.0.2:5001/", new = 2, autoraise=True) # new 2 opens tab while new 1 opens window


########################################################
#                                                      #
#                  Running assembly                    #
#                                                      #
########################################################

# define function for bactflow run in multiprocess for constant realtime steraming.

manager = Manager()
process_status = manager.dict({"pid": None, "running": False, "epoch": 0, "stop_reason": None})
output_queue = Queue()
output_history = manager.list() # to store output history
_current_worker = None
_worker_lock = threading.Lock()
STOP_USER_MSG = "Bactflow has been stopped by the user!"

def _tail_nextflow_log(work_root, output_queue, stop_event, seen_lines):
    """Forward process status/completion lines from .nextflow.log (PTY often misses them)."""
    log_path = os.path.join(work_root, ".nextflow.log")
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
        if os.path.isfile(log_path):
            try:
                with open(log_path, "r", encoding="utf-8", errors="replace") as handle:
                    handle.seek(pos)
                    chunk = handle.read()
                    pos = handle.tell()
            except OSError:
                chunk = ""
            for raw in chunk.splitlines():
                line = clean_log_line(raw)
                if not line or line in seen_lines:
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


def _drain_output_queue():
    """Drop pending stream lines so Stop is not followed by stale exit noise."""
    try:
        while True:
            output_queue.get_nowait()
    except Exception:
        pass


def _bump_epoch():
    """Invalidate any in-flight run_bact worker (Stop or Start takeover)."""
    try:
        epoch = int(process_status.get("epoch") or 0) + 1
    except Exception:
        epoch = 1
    process_status["epoch"] = epoch
    return epoch


def _terminate_worker():
    """Terminate the multiprocessing worker that runs run_bact."""
    global _current_worker
    with _worker_lock:
        worker = _current_worker
        _current_worker = None
    if worker is None:
        return
    try:
        if worker.is_alive():
            worker.terminate()
            worker.join(timeout=1.0)
        if worker.is_alive():
            worker.kill()
    except Exception:
        pass


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


def _stop_pipeline_async(pid):
    """Kill pipeline off the request thread so Stop never blocks the UI."""
    try:
        _kill_process_tree(pid)
        # Final sweep after orphans re-parent.
        _sweep_pipeline_leftovers()
    except Exception as exc:
        print(f"async stop error: {exc}")


def _full_user_stop_async(pid):
    """Stop clicked: kill tree, then terminate worker, then sweep again."""
    try:
        _kill_process_tree(pid)
        _terminate_worker()
        _sweep_pipeline_leftovers()
    except Exception as exc:
        print(f"async user-stop error: {exc}")


def _epoch_is_current(epoch):
    try:
        return int(process_status.get("epoch") or 0) == int(epoch)
    except Exception:
        return False


def run_bact(command, process_status, output_queue, output_history, work_root=None, epoch=0):
    """Function to run bactflow"""

    # Stop/Start may have invalidated this worker before it even started.
    if not _epoch_is_current(epoch):
        return 0

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

    if not _epoch_is_current(epoch):
        stop_event.set()
        if log_thread:
            log_thread.join(timeout=1)
        return 0

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

    if not _epoch_is_current(epoch):
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except Exception:
            try:
                process.kill()
            except Exception:
                pass
        try:
            os.close(master_fd)
        except Exception:
            pass
        stop_event.set()
        if log_thread:
            log_thread.join(timeout=1)
        return 0

    process_status["pid"] = process.pid
    buffer = ""

    while True:
        if not _epoch_is_current(epoch):
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except Exception:
                try:
                    process.kill()
                except Exception:
                    pass
            try:
                process.wait(timeout=1)
            except Exception:
                pass
            break

        still_marked_running = True
        try:
            still_marked_running = bool(process_status["running"])
        except Exception:
            still_marked_running = True

        if (not still_marked_running) and process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except Exception:
                try:
                    process.kill()
                except Exception:
                    pass
            try:
                process.wait(timeout=1)
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

    try:
        os.close(master_fd)
    except Exception:
        pass
    try:
        exit_code = process.wait(timeout=2)
    except Exception:
        try:
            process.kill()
        except Exception:
            pass
        exit_code = process.poll()
        if exit_code is None:
            exit_code = -9

    stop_event.set()
    if log_thread:
        log_thread.join(timeout=2)

    # Decide whether this worker still owns the run lock.
    silence_exit = False
    try:
        if process_status["stop_reason"] == "user":
            silence_exit = True
    except Exception:
        pass
    try:
        if int(process_status["epoch"]) != int(epoch):
            silence_exit = True
    except Exception:
        silence_exit = True
    # SIGKILL / signal exits (e.g. -9) are from Stop — never show them.
    try:
        if exit_code is not None and int(exit_code) < 0:
            silence_exit = True
    except Exception:
        silence_exit = True

    if silence_exit:
        try:
            if int(process_status["epoch"]) == int(epoch):
                process_status["running"] = False
                process_status["pid"] = None
        except Exception:
            pass
        return exit_code

    if exit_code == 0:
        output_queue.put("BactFlow: Nextflow finished successfully.")
    else:
        output_queue.put(f"BactFlow: Nextflow exited with code {exit_code}.")
    process_status["running"] = False
    process_status["pid"] = None
    return exit_code


# run bactflow

@app.route('/run_bactflow', methods=['POST', 'GET'])
def run_bactflow():
    """Startp the BactFlow process"""
    global _current_worker
    if request.method == "POST":
        action = request.args.get("action-assem")
        print("I am here")
        command = None

        if  action == "run":
            # Always take over: Stop may have cleared the lock while an old worker
            # was still about to set running=True / pid — bump epoch so it cannot.
            old_pid = process_status.get("pid")
            epoch = _bump_epoch()
            process_status["stop_reason"] = None
            _clear_run_state()
            _terminate_worker()
            if old_pid:
                threading.Thread(
                    target=_stop_pipeline_async,
                    args=(old_pid,),
                    daemon=True,
                ).start()
            
            
            setup_only = request.form.get("setup_only", 'false')
            cpus = request.form.get('cpus', 1)   
            out_dir = request.form.get('out_dir', './bactflow_out')
            genome_extension = request.form.get('genome_extension', 'fasta')
            checkm_lineag_check = request.form.get('checkm_lineag_check', 'false')
            circle_genome = request.form.get('circle_genome', 'false')  
            tax_class  = request.form.get('tax_class', 'false')
            checkm_db = request.form.get('checkm_db', "")
            gtdbtk_data_path = request.form.get('gtdbtk_data_path', "")
            genome_dir = request.form.get("genome_dir")
            run_checkm = request.form.get("run_checkm", "false")
            run_quast = request.form.get("run_quast", "true")
            resume_run = request.form.get("resume_run", "true")
            run_bakta = request.form.get("bakta_annot", "false")
            bakta_db = request.form.get("bakta_data_path", "")
            command = None

            command = f"""if [ ! -d '{out_dir}' ]; then mkdir -p '{out_dir}'; fi && cd '{base_dir}' && \\
                    nextflow run {base_dir}/main.nf \\
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
                    --genome_dir '{genome_dir}' \\
                    -ansi-log true"""
            if resume_run:
                command = command + " -resume"
            command = with_nextflow_java(command)
                
            output_history[:] = []
            
            back_process = Process(
                target=run_bact,
                args=(command, process_status, output_queue, output_history, base_dir, epoch),
            )
            with _worker_lock:
                _current_worker = back_process
            back_process.start()
        
            command = None
            return "Bactflow started successfully!\n", 200
            
        if action == "help":
            old_pid = process_status.get("pid")
            epoch = _bump_epoch()
            process_status["stop_reason"] = None
            _clear_run_state()
            _terminate_worker()
            if old_pid:
                threading.Thread(
                    target=_stop_pipeline_async,
                    args=(old_pid,),
                    daemon=True,
                ).start()
            command = with_nextflow_java(
                f"cd '{base_dir}' && nextflow run {base_dir}/main.nf --help -ansi-log true"
            )
            output_history[:] = []
            back_process = Process(
                target=run_bact,
                args=(command, process_status, output_queue, output_history, base_dir, epoch),
            )
            with _worker_lock:
                _current_worker = back_process
            back_process.start()
            
            command = None
            return "Bactflow started successfully!\n", 200
        
        if action == "stop":
            pid = process_status.get("pid")
            # Mark user-stop BEFORE kill so the worker never emits "exited with code -9".
            process_status["stop_reason"] = "user"
            _bump_epoch()
            _clear_run_state()
            _drain_output_queue()
            try:
                output_history[:] = []
                output_history.append(STOP_USER_MSG)
            except Exception:
                pass
            threading.Thread(
                target=_full_user_stop_async,
                args=(pid,),
                daemon=True,
            ).start()
            try:
                output_queue.put(STOP_USER_MSG)
            except Exception:
                pass
            return f"{STOP_USER_MSG}\n", 200

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

    if process_status["running"]:
        return {"status": "running"}, 200
    return {"status": "stopped"}, 200


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


@app.route("/sys_stats", methods=["GET"])
def sys_stats():
    try:
        vm = psutil.virtual_memory()
        ncpu = psutil.cpu_count(logical=True) or 1
        cpu = round(psutil.cpu_percent(interval=0.15), 1)
        payload = {
            "cpu_percent": cpu,
            "cpu_count": ncpu,
            "ram_percent": round(vm.percent, 1),
            "ram_used_gb": round(vm.used / (1024 ** 3), 2),
            "ram_total_gb": round(vm.total / (1024 ** 3), 2),
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
        return jsonify(payload)
    except Exception as exc:
        return jsonify({
            "cpu_percent": 0,
            "cpu_count": 1,
            "ram_percent": 0,
            "ram_used_gb": 0,
            "ram_total_gb": 0,
            "running": False,
            "error": str(exc),
        })

#Stream    
@app.route('/stream_bactflow', methods = ['POST', 'GET'])
def stream_bactflow():

  
    def generate():
        # fitst show the history
        for line in list(output_history):
            yield f"data: {line}\n\n"

        # now stream new output
        idle_rounds = 0
        while True:
            running = False
            try:
                running = bool(process_status["running"])
            except Exception:
                running = False

            try:
                line = output_queue.get(timeout=0.4)
                idle_rounds = 0
                yield f"data: {line.strip()}\n\n"
            except Exception:
                idle_rounds += 1
                if (not running) and idle_rounds >= 2:
                    break
                continue

        try:
            if process_status.get("stop_reason") == "user":
                yield f"data: {STOP_USER_MSG}\n\n"
            else:
                yield "data: Process completed\n\n"
        except Exception:
            yield "data: Process completed\n\n"

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

# quast
@app.route("/check-quast", methods = ["POST"])  
def check_quast():
    out_dir = request.form.get("out_dir")
    if not out_dir:
        return jsonify({"exists":False, "error": "Missing out_dir"}), 400
    quast_path = os.path.join(out_dir, "quast_stat")
    if os.path.exists(quast_path):
        return jsonify({"exists": True})
    else:
        return jsonify({"exists": False})


@app.route("/quast-report", methods = ["POST"])
def quast_report():
    out_dir = request.form.get("out_dir")
    quast_path = os.path.join(out_dir, "quast_stat")

    if os.path.exists(quast_path):
        return send_from_directory(quast_path, "report.html")
    else:
        return jsonify({"error": "QUAST report does not exist"}), 404

@app.route("/contig-report", methods = ["POST"])
def contig_report():
    out_dir = request.form.get("out_dir")
    quast_path = os.path.join(out_dir, "quast_stat")
    if os.path.exists(quast_path):
        return send_from_directory(quast_path,"icarus_viewers/contig_size_viewer.html")
    else:
        return jsonify({"error": "Contig size report does not exist"}), 404

# Bakta
@app.route("/check-bakta", methods = ["POST"])  
def check_bakta():
    out_dir = request.form.get("out_dir")
    gt = request.form.getlist("gene_type")
    gene_type = ",".join(gt) 
    
    if not out_dir:
        return jsonify({"exists":False, "error": "Missing out_dir"}), 400
    bakta_path = os.path.join(out_dir, "bakta_out")
    if os.path.exists(bakta_path):
        try:

            command = f"""
            cp {bakta_path}/*/*.tsv {bakta_path}
            rm {bakta_path}/*inference.tsv {bakta_path}/*hypotheticals.tsv 
            mkdir -p {bakta_path}/genes && mv {bakta_path}/*.tsv {bakta_path}/genes
           
            python {base_dir}/gene_counter_bakta.py -d {bakta_path}/genes -t {gene_type} -o {out_dir}/gene_count
            """
            
            gene_count = os.path.join(out_dir, "gene_count.tsv")
            subprocess.run(command, shell=True, text=True, check=True, capture_output=True)
            if not os.path.exists(gene_count):
                subprocess.run(command, shell=True, text=True, check=True, capture_output=True)

            df = pd.read_csv(gene_count, sep="\t")

            
            if os.path.exists(gene_count) and len(gene_count) != 0:
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
                # "message": result.stdout.strip(), #remove whitespace 
                # "error": result.stderr.strip(), 
                "count_tab": render_template_string(
                    table_html, id="bakta-tab", tabnumber = "Table 1: Count table of annotated genes by Bakta.", tabcaption="This is the abundance of all genes (selected based on gene type input) at all genomes.", table=table_data
                )
                })
        except subprocess.CalledProcessError as e:
            return jsonify({
                "exists": True,
                "error": f"Gene counter failed with error code {e.returncode}",
                "stderr": e.stderr.strip(),
                "stdout": e.stdout.strip()

                }), 500

            
    else:
        return jsonify({"exists": False, "error": "Bakta annotation directory doesn't exist!"}), 400

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


    
@app.route("/taxa-report", methods = ["POST"])
def taxa_report():
    out_dir = request.form.get("out_dir")
    tax_file = os.path.join(out_dir, "gtdbtk_out/classify/gtdbtk.bac120.summary.tsv")

    if os.path.exists(tax_file):
        df = pd.read_csv(tax_file, sep = "\t")
        if all(col in df.columns for col in ['user_genome', 'classification', 'closest_genome_ani']):
            df  = df[['user_genome', 'classification', 'closest_genome_ani']]

            data = df.to_dict(orient="records")

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
                 "taxa_table": render_template_string(
                    table_html, id="taxa-tab", tabnumber = "Table 2: Taxonomy classification table of genomes.", tabcaption="This classification was done by GDTBdk", table=data
                )
                })
        else:
            return jsonify({"exists": False, "error": "Missing expected columns in taxonomy.tsv"}), 400
    else:
        return jsonify({"exists": False, "error": "taxonomy.tsv file not found"}), 404

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

    snp_files = glob.glob(os.path.join(out_dir, "snps", "**", "*.vcf"), recursive=True)
    if snp_files:
        return jsonify({"exists": True})
    return jsonify({"exists": False, "error": "SNP finder finished but no VCF files were found."})

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

    vcf_files = glob.glob(os.path.join(out_dir, "vcs", "**", "*.vcf"), recursive=True)
    if vcf_files:
        return jsonify({"exists": True})
    return jsonify({"exists": False, "error": "Variant calling finished but no VCF files were found."})

# abundance
@app.route("/abund-run", methods = ["POST"])
def abund_finder():
    out_dir = request.form.get("out_dir")
    gene_files = request.form.get("gene_files")
    prevalance = "false"
    gene_type = ",".join(_selected_gene_types())
    width = request.form.get("plot_width") or "10"
    height = request.form.get("plot_height") or "10"

    
    enzyme_file = request.form.get("enzyme_loc")
    if not out_dir:
        return jsonify({"exists": False, "error": "Output directory is missing."}), 400
    if not gene_files:
        return jsonify({"exists": False, "error": "Gene annotation directory is missing."}), 400
    if not enzyme_file:
        return jsonify({"exists": False, "error": "Enzyme file path is missing."}), 400
    output = os.path.join(out_dir, "strain_finder")
    count_tab = os.path.join(output, "abundance.tsv")
    plot = os.path.join(output, "requested_genes_abundance.png")
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

    if os.path.exists(count_tab) and os.path.exists(plot):
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
        
        with open(plot, "rb") as img_file:
            img_base = base64.b64encode(img_file.read()).decode("utf-8")
      
        return jsonify({
                "exists": True, 
                 "abund_table": render_template_string(
                    table_html, id="abund-tab", tabnumber = "Table 3: Abundance table of requested genes.", table=data
                ),
                "plot_abund" : f"data:image/png;base64,{img_base}"
                })
    else:
        return jsonify({"exists": False, "error": "Missing expected table and plot"}), 400



# prevalance
@app.route("/prev-run", methods = ["POST"])
def prev_finder():
    out_dir = request.form.get("out_dir")
    gene_files = request.form.get("gene_files")
    prevalance = "true"
    gene_type = ",".join(_selected_gene_types())
    
    enzyme_file = request.form.get("enzyme_loc")
    if not out_dir:
        return jsonify({"exists": False, "error": "Output directory is missing."}), 400
    if not gene_files:
        return jsonify({"exists": False, "error": "Gene annotation directory is missing."}), 400
    if not enzyme_file:
        return jsonify({"exists": False, "error": "Enzyme file path is missing."}), 400
    output = os.path.join(out_dir, "strain_finder")
    count_tab = os.path.join(output, "prevalance.tsv")
    plot = os.path.join(output, "requested_genes_prevalence.png")
    width = request.form.get("plot_width") or "10"
    height = request.form.get("plot_height") or "10"

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

    if os.path.exists(count_tab) and os.path.exists(plot):
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
        
        with open(plot, "rb") as img_file:
            img_base = base64.b64encode(img_file.read()).decode("utf-8")
      
        return jsonify({
                "exists": True, 
                 "prev_table": render_template_string(
                    table_html, id="prev-tab", tabnumber = "Table 4: Prevalance table of requested genes.", table=data
                ),
                "plot_prev" : f"data:image/png;base64,{img_base}"
                })
    else:
        return jsonify({"exists": False, "error": "Missing expected table and plot"}), 400




if __name__ == '__main__':
    Timer(1, open_browser).start()
    app.run(debug=True, port=5001, host="0.0.0.0", use_reloader=False, threaded=True)#set use_reloader to true during developement 
