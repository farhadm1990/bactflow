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


def clean_log_line(line):
    line = ANSI_RE.sub("", line or "")
    return line.replace("\r", "").strip()





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
        print("I am here")
        command = None

        if  action == "run":
            if process_status["running"]:
                return "Bactflow is already running!", 400
            
            
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
                args=(command, process_status, output_queue, output_history, base_dir),
            )
            back_process.start()
        
            command = None
            return "Bactflow started successfully!\n", 200
            
        if action == "help":
            command = with_nextflow_java(
                f"cd '{base_dir}' && nextflow run {base_dir}/main.nf --help -ansi-log true"
            )
            output_history[:] = []
            back_process = Process(
                target=run_bact,
                args=(command, process_status, output_queue, output_history, base_dir),
            )
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
        for line in output_history:
            yield f"data: {line}\n\n"

        # now stream new output
        while process_status['running']  or not output_queue.empty():
            try:
                line = output_queue.get(timeout=0.5)# wait for output
                yield f"data: {line.strip()}\n\n"
            except Exception:
                if not process_status['running']:
                    break
                continue

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

@app.route("/check-bakta-ready", methods=["POST"])
def check_bakta_ready():
    """Lightweight check for circular-plot input (.gbk files in bakta_out)."""
    out_dir = request.form.get("out_dir", "").strip()
    if not out_dir:
        return jsonify({
            "ready": False,
            "message": "Set an output directory first."
        })

    bakta_path = os.path.join(out_dir, "bakta_out")
    if not os.path.isdir(bakta_path):
        return jsonify({
            "ready": False,
            "message": "No Bakta output yet. Run BactFlow with Bakta enabled."
        })

    gbk_files = glob.glob(os.path.join(bakta_path, "**", "*.gbk"), recursive=True)
    if gbk_files:
        return jsonify({"ready": True, "gbk_count": len(gbk_files)})

    return jsonify({
        "ready": False,
        "message": "Bakta output exists but no .gbk files were found. Check the Bakta run log for errors."
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
        "figsize": int(request.form.get("figsize", 200)),
        "interval": int(request.form.get("interval", 3)),
        "f_color": request.form.get("f_color", "#1E90FF"),
        "r_color": request.form.get("r_color", "#FF7261"),
    }

    last_params = {}


    if os.path.exists(params_file):
        try:
            with open(params_file, 'r') as f:
                if os.path.getsize(params_file) > 0:
                    last_params = json.load(f)

        except (json.JSONDecodeError, IOError) as e:
            print(f"Error reading params file for the plot: {e}")
    
    needs_regen = not os.path.exists(crc_plt) or params != last_params

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
        bash -c "source $(conda info --base)/etc/profile.d/conda.sh && \
        conda activate bactflow && \
        python3 {base_dir}/circular_plotter.py -d {gbk_dir} -o {out_dir} \
        {'--add_gc' if params['add_gc'] == 'True' else ''} \
        {'--add_skew' if params['add_skew'] == 'True' else ''} \
        --dpi {params['dpi']} \
        --interval {params['interval']} \
        --f_color '{params['f_color']}' \
        --r_color '{params['r_color']}'"
        """

        try:
            subprocess.run(
                command, shell=True, check=True, text=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            print(" am running the command ")
        except Exception as e:
            return jsonify({"plot": False, "error": str(e)}), 200

    if os.path.exists(crc_plt):
        with open(crc_plt, "rb") as img_file:
            img_base = base64.b64encode(img_file.read()).decode("utf-8")
        return jsonify({"plot": f"data:image/png;base64,{img_base}"})
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
    print(reference)
    command = f"""
    conda activate bactflow 
    {base_dir}/variant_finder.sh  -r {reference} -g {genomes} -o {out_dir}/snps -t {cpus}
    """

    subprocess.run(command, shell=True, check=True, text=True)

    snp_file = glob.glob(os.path.join(out_dir, "snp/*/*.vcf"))

    if os.path.exists(snp_file):
        return jsonify({"exists": True})
    else:
        return jsonify({"exists": False})

@app.route("/svs-finder", methods = ["POST"])
def svs_finder():
    out_dir = request.form.get("out_dir")
    genomes = request.form.get("genome_dir")
    reference = request.form.get("ref_genome")
    cpus = request.form.get("cpus")
  
    command = f"""
    conda activate bactflow 
    {base_dir}/vc_medaka.sh  -r {reference} -g {genomes} -o {out_dir}/vcs -c {cpus}
    """

    subprocess.run(command, shell=True, check=True, text=True)

    snp_file = glob.glob(os.path.join(out_dir, "vcs/*/*.vcf"))

    if os.path.exists(snp_file):
        return jsonify({"exists": True})
    else:
        return jsonify({"exists": False})

# abundance
@app.route("/abund-run", methods = ["POST"])
def abund_finder():
    out_dir = request.form.get("out_dir")
    gene_files = request.form.get("gene_files")
    prevalance = "false"
    gene_type = request.form.get("gene_type").lower()
    width = request.form.get("plot_width")
    height = request.form.get('plot_height')

    
    output = os.path.join(out_dir, "strain_finder")
    enzyme_file = request.form.get("enzyme_loc")
    count_tab = os.path.join(output, "abundance.tsv")
    plot = os.path.join(output, "requested_genes_abundance.png")
    print(f"this is widht {width} and heigth {height}")
    command = f"""
    {base_dir}/strain_finder.sh -d {gene_files} -g {gene_type} -f abundance -o {output} -e {enzyme_file}  -c false -p {prevalance} -w {width} -l {height}

    """
    
    
    subprocess.run(command, shell=True, text=True, check=True)

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
    gene_type = request.form.get("gene_type").lower()
    
    output = os.path.join(out_dir, "strain_finder")
    enzyme_file = request.form.get("enzyme_loc")
    count_tab = os.path.join(output, "prevalance.tsv")
    plot = os.path.join(output, "requested_genes_prevalence.png")
    width = request.form.get("plot_width")
    height = request.form.get('plot_height')

    command = f"""
    
    {base_dir}/strain_finder.sh -d {gene_files} -g {gene_type} -f prevalance -o {output} -e {enzyme_file}  -c false -p {prevalance} -w {width} -l {height}

    """
   
    subprocess.run(command, shell=True, text=True, check=True)

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
    app.run(debug = True, port = 5001, host = "0.0.0.0",  use_reloader = False)#set use_reloader to true during developement 
