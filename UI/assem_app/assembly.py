#!/usr/bin/env python3

import sys, subprocess, importlib, os

from flask import Flask, send_from_directory, render_template, request, redirect, Response, stream_with_context, jsonify

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
import threading
from threading import Timer
import webbrowser
import re


base_dir = os.path.abspath(os.path.dirname(__file__))# we can have access to all files from everywhere
app = Flask(__name__, 
            template_folder = os.path.join(base_dir, "templates"),
            static_folder = os.path.join(base_dir, "static"))

# Nextflow accepts Java 17-24. The bactflow conda env currently ships OpenJDK 25.
NF_JAVA_SETUP = r"""
unset JAVA_CMD
java_major() {
  "$1" -version 2>&1 | sed -n 's/.*version "\([0-9]\+\).*/\1/p' | head -n1
}
pick_nf_java() {
  local candidate major java_bin
  for candidate in \
    ${JAVA_HOME:+"$JAVA_HOME/bin/java"} \
    "$(command -v java 2>/dev/null)" \
    "$HOME/.sdkman/candidates/java/current/bin/java" \
    /usr/lib/jvm/java-21-openjdk-amd64/bin/java \
    /usr/lib/jvm/java-17-openjdk-amd64/bin/java
  do
    [ -n "$candidate" ] && [ -x "$candidate" ] || continue
    major=$(java_major "$candidate")
    if [ -n "$major" ] && [ "$major" -ge 17 ] && [ "$major" -le 24 ]; then
      java_bin=$(readlink -f "$candidate" 2>/dev/null || echo "$candidate")
      export JAVA_CMD="$java_bin"
      export JAVA_HOME="$(dirname "$(dirname "$java_bin")")"
      export PATH="$(dirname "$java_bin"):$PATH"
      return 0
    fi
  done
  return 1
}
if ! pick_nf_java; then
  echo "ERROR: Nextflow needs Java 17-24. The bactflow env Java is too new (OpenJDK 25)." >&2
  exit 1
fi
echo "Using Java for Nextflow: $JAVA_CMD"
"$JAVA_CMD" -version
export NXF_ANSI_LOG=false
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

def run_bact(command, process_status, output_queue, output_history):
    """Function to run bactflow"""


    process_status["running"]=True
    process = subprocess.Popen(
        command, 
        shell=True,
        stdout=subprocess.PIPE, 
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1
    )
    process_status["pid"]=process.pid
    output_history[:] = []
    history_cap = 1500
    for line in iter(process.stdout.readline, ""):
        line = clean_log_line(line)
        if not line:
            continue
        output_queue.put(line)
        output_history.append(line)
        if len(output_history) >= history_cap + 250:
            del output_history[0:250]

    process.stdout.close()
    process.wait()
    process_status["running"] = False
    process_status['pid'] = None


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
            pacbio_read_type = request.form.get('pacbio_read_type', 'hifi')
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
             
            command = f"""if [ ! -d '{out_dir}' ]; then mkdir -p '{out_dir}'; fi && 
                    nextflow run {base_dir}/main.nf \
                    --setup_only {setup_only} \
                    --fastq_dir '{fastq_dir}' \
                    --concat_reads {concat_reads} \
                    --extension {extension}\
                    --cpus {cpus} \
                    --coverage_filter {coverage_filter} \
                    --coverage {coverage} \
                    --genome_size {genome_size} \
                    --out_dir '{out_dir}' \
                    --tensor_batch {tensor_batch} \
                    --nanofilter {nanofilter} \
                    --min_length {min_length} \
                    --min_quality {min_quality} \
                    --medaka_polish {medaka_polish} \
                    --basecaller_model {basecaller_model} \
                    --genome_extension {genome_extension} \
                    --checkm_lineag_check {checkm_lineag_check} \
                    --run_flye {run_flye} \
                    --circle_genome {circle_genome} \
                    --run_unicycler {run_unicycler} \
                    --run_spades {run_spades} \
                    --run_pacbio {run_pacbio} \
                    --short_read_dir '{short_read_dir}' \
                    --pacbio_read_type {pacbio_read_type} \
                    --tax_class {tax_class} \
                    --bakta_annot {bakta_annot} \
                    --bakta_db {bakta_db} \
                    --run_checkm {run_checkm} \
                    --checkm_db {checkm_db} \
                    --gtdbtk_data_path {gtdbtk_data_path} \
                    --run_quast {run_quast} \
                    --genome_dir {genome_dir} \
                    -ansi-log false"""
            if resume_run:
                command = command + " -resume"
            else:
                command = command
            command = with_nextflow_java(command)
                
            output_history[:] = []
            
            back_process = Process(target=run_bact, args=(command, process_status, output_queue, output_history))
            back_process.start()
        
            command = None
            return "Bactflow started successfully!\n", 200
            
        if action == "help":
            command = with_nextflow_java(f"nextflow run {base_dir}/main.nf --help -ansi-log false")
            output_history[:] = []
            back_process = Process(target=run_bact, args=(command, process_status, output_queue, output_history))
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
# prog = {"completed": 0}

# @app.route("/progress", methods = ["POST", "GET"])
# def progress():
    
#     if request.method == "POST":
#         out_dir = request.form.get('out_dir', './bactflow_out')
#         genome_dir = request.form.get('genome_dir')


#         pars =  None

#         if not genome_dir:
#             pars = {
#             "run_flye": True if request.form.get("run_flye") == "true" else False,
#             "run_unicycler": True if request.form.get("run_unicycler") == "true" else False,
#             "run_spades": True if request.form.get("run_spades") == "true" else False,
#             "run_megahit": True if request.form.get("run_megahit") == "true" else False,
#             "run_quast": True if request.form.get("run_quast") == "true" else False,
#             "circle_genome": True if request.form.get("circle_genome") == "true" else False,
#             "tax_class": True if request.form.get("tax_class") == "true" else False,
#             "prok_annot": True if request.form.get("prok_annot") == "true" else False,
#             "run_checkm": True if request.form.get("run_checkm") == "true" else False,
#             "run_bakta": True if request.form.get("run_bakta") == "true" else False,
#             }
#         else:
#             pars = {
#             "run_quast": True if request.form.get("run_quast") == "true" else False,
#             "tax_class": True if request.form.get("tax_class") == "true" else False,
#             "prok_annot": True if request.form.get("prok_annot") == "true" else False,
#             "run_checkm": True if request.form.get("run_checkm") == "true" else False,
#             "run_bakta": True if request.form.get("run_bakta") == "true" else False,
#             }

       
#         def out_dir_watcher(out_dir):
            
#             global prog
            
#             expected_counts = sum(value is True for value in pars.values())

#             pr_stat = process_status["running"]

#             if expected_counts > 0:
#                 while pr_stat:
#                     try:
#                         current_count = len(os.listdir(out_dir)) if os.path.exists(out_dir) else 0
#                         prog = 0
#                         prog = {"completed": f"{current_count/expected_counts*100}"}

#                         if current_count >= expected_counts:
#                             pr_stat = False
#                     except Exception as e:
#                         print(f"Error watching the output directory: {e}")
#                     time.sleep(2)
#         wather_thread = threading.Thread(target=out_dir_watcher, args=(out_dir,), daemon=True)
#         wather_thread.start() 

#     return jsonify(prog) 

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
    quast_check = True if os.path.exists(f"{out_dir}/quast_stat") else False
    if quast_check:
        return send_from_directory(f"{out_dir}/quast_stat","report.html")

@app.route("/contig-report", methods = ["POST"])
def contig_report():
    out_dir = request.form.get("out_dir")
    quast_check = True if os.path.exists(f"{out_dir}/quast_stat") else False
    if quast_check:
        return send_from_directory(f"{out_dir}/quast_stat","icarus_viewers/contig_size_viewer.html")

if __name__ == '__main__':
    Timer(1, open_browser).start()
    app.run(debug = True, port = 5002, host = "0.0.0.0",  use_reloader = False)#set use_reloader to true during developement 
