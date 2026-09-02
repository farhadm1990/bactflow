#!/usr/bin/env python3

import sys, subprocess, importlib, os
import pty
import select

from flask import Flask, send_from_directory, render_template, request, redirect, Response, stream_with_context, jsonify
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
import threading
from threading import Timer
import webbrowser
import re


base_dir = os.path.abspath(os.path.dirname(__file__))# we can have access to all files from everywhere
IN_DOCKER = os.environ.get("BACTFLOW_IN_DOCKER") == "1" or os.path.exists("/.dockerenv")
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


def clean_log_line(line):
    line = ANSI_RE.sub("", line or "")
    return line.replace("\r", "").strip()





#assembly
@app.route('/',  methods =['GET', 'POST'])
def assembly():

    return render_template('assembly.html')


def open_browser():
    webbrowser.open(url="http://127.0.0.1:5002/", new = 2, autoraise=True) # new 2 opens tab while new 1 opens window


########################################################
#                                                      #
#                  Running assembly                    #
#                                                      #
########################################################

# define function for bactflow run in multiprocess for constant realtime steraming.

manager = Manager()
process_status = manager.dict({"pid": None, "running": False})
output_queue = Queue()
output_history = manager.list() # to store output history 

def _tail_nextflow_log(work_root, output_queue, stop_event, seen_lines):
    """Forward completion lines from .nextflow.log (PTY stream often misses them)."""
    log_path = os.path.join(work_root, ".nextflow.log")
    pos = 0
    skip_re = re.compile(
        r"\[Task monitor\]|TaskPollingMonitor|\bDEBUG\b",
        re.IGNORECASE,
    )
    done_re = re.compile(
        r"✔|Executor finished|\[100%\].*of",
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
                if done_re.search(line) or (
                    re.search(r"\[[0-9a-f/]+\]", line, re.I)
                    and re.search(r"process\s+>", line, re.I)
                    and re.search(r"\[\s*\d+%\]", line)
                ):
                    seen_lines.add(line)
                    output_queue.put(line)
        stop_event.wait(0.8)


def run_bact(command, process_status, output_queue, output_history, work_root=None):
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
        )
    finally:
        os.close(slave_fd)

    process_status["pid"] = process.pid
    buffer = ""

    while True:
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
    if request.method == "POST":
        action = request.args.get("action-assem")
        
        command = None

        if  action == "run":
            if process_status["running"]:
                return "Bactflow is already running!", 400
            
            
            setup_only = request.form.get("setup_only", 'false')
            fastq_dir = request.form.get('fastq_dir', './')

            # fastq_file = request.files.get('fastq_file')
            # if fastq_file:
            #     save_path = os.path.join("/tmp", fastq_file.filename)
            #     fastq_dir.save(save_path)
            concat_reads = request.form.get("concat_reads", "true")
            extension = request.form.get('extension', '.fastq.gz')
            cpus = request.form.get('cpus', 1)   
            coverage_filter = request.form.get('coverage_filter', 'false')
            coverage = request.form.get('coverage', 40)
            genome_size = request.form.get('genome_size', 6)
            out_dir = request.form.get('out_dir', './bactflow_out')
            tensor_batch = request.form.get('tensor_batch', 200)
            nanofilter = request.form.get('nanofilter', 'true')
            min_length = request.form.get('min_length', 1000)
            min_quality = request.form.get('min_quality', 16)
            medaka_polish = request.form.get('medaka_polish', 'false')
            basecaller_model = request.form.get('basecaller_model', 'r1041_e82_400bps_hac_v4.2.0')
            genome_extension = request.form.get('genome_extension', 'fasta')
            checkm_lineag_check = request.form.get('checkm_lineag_check', 'false')
            run_flye = request.form.get('run_flye', 'false')
            circle_genome = request.form.get('circle_genome', 'false')  
            run_unicycler = request.form.get('run_unicycler', 'false')
            run_spades = request.form.get('run_spades', 'false')
            run_pacbio = request.form.get('run_pacbio', 'false')
            short_read_dir = request.form.get('short_read_dir', '')
            ont_read_type = request.form.get('ont_read_type', 'nano-raw')
            pacbio_read_type = request.form.get('pacbio_read_type', 'pacbio-hifi')
            tax_class  = request.form.get('tax_class', 'false')
            checkm_db = request.form.get('checkm_db', "")
            gtdbtk_data_path = request.form.get('gtdbtk_data_path', "")
            genome_dir = request.form.get("genome_dir")
            run_checkm = request.form.get("run_checkm", "false")
            run_quast = request.form.get("run_quast", "true")
            resume_run = request.form.get("resume_run", "true")
            bakta_annot = request.form.get("bakta_annot", "false")
            bakta_db = request.form.get("bakta_data_path", "")
            command = None
            print(f"this is bakta_annot {bakta_annot}")

            if fastq_dir:
                fastq_dir = os.path.abspath(fastq_dir)
            if short_read_dir:
                short_read_dir = os.path.abspath(short_read_dir)
            if out_dir:
                out_dir = os.path.abspath(out_dir)

            if str(run_unicycler).lower() == "true" and not str(short_read_dir).strip():
                return "Unicycler hybrid assembly requires a short-read FASTQ directory.\n", 400
            if str(run_unicycler).lower() == "true" and not str(fastq_dir).strip():
                return "Unicycler hybrid assembly requires a long-read FASTQ directory.\n", 400
            if str(run_spades).lower() == "true" and not str(fastq_dir).strip():
                return "SPAdes requires an Illumina paired-end FASTQ directory.\n", 400
             
            command = f"""if [ ! -d '{out_dir}' ]; then mkdir -p '{out_dir}'; fi && cd '{base_dir}' && \\
                    nextflow run {base_dir}/main.nf \\
                    --setup_only {setup_only} \\
                    --fastq_dir '{fastq_dir}' \\
                    --concat_reads {concat_reads} \\
                    --extension {extension}\\
                    --cpus {cpus} \\
                    --coverage_filter {coverage_filter} \\
                    --coverage {coverage} \\
                    --genome_size {genome_size} \\
                    --out_dir '{out_dir}' \\
                    --tensor_batch {tensor_batch} \\
                    --nanofilter {nanofilter} \\
                    --min_length {min_length} \\
                    --min_quality {min_quality} \\
                    --medaka_polish {medaka_polish} \\
                    --basecaller_model {basecaller_model} \\
                    --genome_extension {genome_extension} \\
                    --checkm_lineag_check {checkm_lineag_check} \\
                    --run_flye {run_flye} \\
                    --ont_read_type {ont_read_type} \\
                    --circle_genome {circle_genome} \\
                    --run_unicycler {run_unicycler} \\
                    --run_spades {run_spades} \\
                    --run_pacbio {run_pacbio} \\
                    --short_read_dir '{short_read_dir}' \\
                    --pacbio_read_type {pacbio_read_type} \\
                    --tax_class {tax_class} \\
                    --bakta_annot {bakta_annot} \\
                    --bakta_db {bakta_db} \\
                    --run_checkm {run_checkm} \\
                    --checkm_db {checkm_db} \\
                    --gtdbtk_data_path {gtdbtk_data_path} \\
                    --run_quast {run_quast} \\
                    --genome_dir {genome_dir} \\
                    -ansi-log true"""
            if resume_run:
                command = command + " -resume"
            command = command + f"""
                    NF_EXIT=$?
                    if [ $NF_EXIT -eq 0 ]; then
                      echo "BactFlow: cleaning temporary Nextflow work files..."
                      rm -rf '{out_dir}/.nextflow-work' '{base_dir}/work' 2>/dev/null || true
                      "$BACTFLOW_NEXTFLOW_BIN" clean -f 2>/dev/null || true
                    fi
                    exit $NF_EXIT"""
            command = with_nextflow_java(command)
                
            output_history[:] = []
            
            back_process = Process(target=run_bact, args=(command, process_status, output_queue, output_history, base_dir))
            back_process.start()
        
            command = None
            return "Bactflow started successfully!\n", 200
            
        if action == "help":
            command = with_nextflow_java(f"cd '{base_dir}' && nextflow run {base_dir}/main.nf --help -ansi-log true")
            output_history[:] = []
            back_process = Process(target=run_bact, args=(command, process_status, output_queue, output_history, base_dir))
            back_process.start()
            
            command = None
            return "Bactflow started successfully!\n", 200
        
        if action == "stop":
            if process_status["running"] and process_status["pid"]:
                try:
                    parent = psutil.Process(process_status['pid'])
                    for child in parent.children(recursive=True):#terminate child first
                        child.terminate()
                    parent.terminate()
                    gone, still_alive = psutil.wait_procs([parent], timeout=5)
                    for p in still_alive:
                        p.kill()
                    
                

                    # os.kill(process_status["pid"], signal.SIGTERM) 
                    process_status["running"] = False
                    process_status["pid"] = None
                    
                    return "Bactflow stopped successfully!\n", 200
                except psutil.NoSuchProcess:
                    
                    return "Process not found!\n", 400
                except Exception as e:
                    return f"Error stopping process: {str(e)}\n", 500
            else:
                return "No running process to stop.\n", 400

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


def _build_sys_stats():
    vm = psutil.virtual_memory()
    ncpu = psutil.cpu_count(logical=True) or 1
    payload = {
        "cpu_percent": round(psutil.cpu_percent(interval=None), 1),
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
    return payload


_SYS_STATS = {}
_SYS_STATS_LOCK = threading.Lock()
_SYS_STATS_JSON = os.path.join(base_dir, "static", "sys_stats.json")


def _write_sys_stats_file(payload):
    try:
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

        yield "data: Process completed\n\n"
 

    return Response(generate(), content_type='text/event-stream')

def find_quast_dir(out_dir):
    """Use the single combined QUAST report in quast_stat."""
    if not out_dir:
        return None
    path = os.path.join(out_dir, "quast_stat")
    if os.path.isdir(path) and os.path.isfile(os.path.join(path, "report.html")):
        return path
    return None

@app.route("/check-quast", methods = ["POST"])  
def check_quast():
    out_dir = request.form.get("out_dir")
    if not out_dir:
        return jsonify({"exists":False, "error": "Missing out_dir"}), 400
    quast_path = find_quast_dir(out_dir)
    if quast_path:
        return jsonify({"exists": True, "quast_dir": os.path.basename(quast_path)})
    return jsonify({"exists": False})


@app.route("/quast-report", methods = ["POST"])
def quast_report():
    out_dir = request.form.get("out_dir")
    quast_path = find_quast_dir(out_dir)
    if quast_path:
        return send_from_directory(quast_path, "report.html")
    return jsonify({"error": "QUAST report not found"}), 404

@app.route("/contig-report", methods = ["POST"])
def contig_report():
    out_dir = request.form.get("out_dir")
    quast_path = find_quast_dir(out_dir)
    if quast_path:
        return send_from_directory(quast_path, "icarus_viewers/contig_size_viewer.html")
    return jsonify({"error": "QUAST contig viewer not found"}), 404

if __name__ == '__main__':
    skip_browser = IN_DOCKER or os.environ.get("BACTFLOW_NO_BROWSER") == "1"
    if not skip_browser:
        Timer(1, open_browser).start()
    app.run(
        debug=not skip_browser,
        port=5002,
        host="0.0.0.0",
        use_reloader=False,
        threaded=True,
    ) 
