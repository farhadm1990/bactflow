
# installing pkgs

import sys, subprocess, importlib, os

def install_pks(pks):
    """installing pkgs"""
 
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install"] + pks)
    except subprocess.CalledProcessError as e:
        print(f"Failed to install packges: {pks}")
        sys.exit(1)

required_packages = {
    'flask': 'Flask',
    'plotly': 'plotly',
    'biopython': 'Bio',
    'flask_sqlalchemy': 'Flask-SQLAlchemy',
    'flask_migrate': 'Flask-Migrate',
    'pandas': 'pandas',
    'numpy': 'numpy',
    'psutil': 'psutil'
}   

missing_pkg = []

for module_name, package_name in required_packages.items():
    try:
        importlib.import_module(module_name)
    except ImportError:
        print(f"Package {package_name} not found. It will be installed now...!")
        missing_pkg.append(package_name)

if missing_pkg:
    install_pks(missing_pkg)
    print("Restarting the script to apply changes...")
    # os.execv(sys.executable, [sys.executable] + sys.argv)


from flask import Flask, render_template, request, redirect, Response, stream_with_context, jsonify, render_template_string
import plotly.express as px
import concurrent.futures
import plotly as py
from Bio import SeqIO
from flask_sqlalchemy import SQLAlchemy
import json
import gzip
from datetime import datetime, timezone
import smtplib # for emails 
import pandas as pd
import numpy as np
import os
# from dotenv import load_dotenv
from flask_migrate import Migrate
import subprocess
import sys
from multiprocessing import Process, Manager, Queue
import time
import signal
import psutil  # for better process control
from threading import Timer
import webbrowser


base_dir = os.path.abspath(os.path.dirname(__file__))
app = Flask(__name__, 
            template_folder = os.path.join(base_dir, "templates"),
            static_folder = os.path.join(base_dir, 'static'))

manager = Manager()
process_status = manager.dict({"running": False})
output_queue = Queue()

def run_bact(command, process_status, output_queue):
    """To run bactflow in mutliprocess"""
    process_status["running"] = True
    process = subprocess.Popen(
        command, 
        shell=True, 
        stdout= subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )

    process_status["pid"] = process.pid
    for line in iter(process.stdout.readline, ""):
        output_queue.put(line.strip())
        time.sleep(0.1)

    process.stdout.close()
    process.wait()
    process_status["running"] = False



#pre assembly
@app.route('/',  methods =['GET', 'POST'])
def preassembly():
    
    return render_template('pre-assembly.html')

def open_browser():
    webbrowser.open(url="http://127.0.0.1:5000/", new=2)

# First check if bactflow is installed

@app.route('/bactflow-check', methods = ['GET'])
def bact_check():
    cmd = f""" 
    
    source $(conda info --base)/etc/profile.d/conda.sh
    conda config --set report_errors false

    if command -v mamba &>/dev/null; then
        env_list=$(mamba env list --json)
    else
        env_list=$(conda env list --json)
    fi

    echo "$env_list"
    """
    try:

        result = subprocess.run(cmd, shell=True, text = True, executable="/bin/bash", capture_output=True)

        output = result.stdout.strip()
        error = result.stderr.strip()

        if error:
            return jsonify({"status": "error", "message": error}), 500
        
        envs_data = json.loads(output)
        envs = envs_data.get("envs", [])
        is_installed = any("bactflow" in env for env in envs)
        
        
        return jsonify({"status": "success", "installed": is_installed}), 200
        
    except Exception as e:
        return jsonify({"status": "error", "message": str(e)}), 500
 


# install bactflow

@app.route('/install-bactflow', methods = ['POST', 'GET'])
def install_bactflow():
    if request.method == 'POST':

        out_dir = request.form.get("out_dir")
        if process_status["running"]:
            return jsonify({"message": "Installation is already in progress", "running": True}), 400
        process_status["running"] = True
            


        command = f"""
        
        nextflow run {base_dir}/main.nf --setup_only true --out_dir {out_dir}
        """
    
        back_process = Process(target=run_bact, args=(command, process_status, output_queue))
        back_process.start()
      
        return jsonify({"message": "Installation started", "running": True}), 200
    elif request.method == 'GET':
        return jsonify({"running": process_status["running"]})

# concat reads and read stat route
@app.route('/ls-fastq', methods = ['POST', 'GET'])
def ls_fastq():
    if request.method == 'POST':
        fastq_dir = request.form.get("fastq_dir")
        extension = request.form.get("extension")
        out_dir = request.form.get("out_dir")
        cpus = request.form.get("cpus")
        concater = request.form.get("concat_reads")
   
        if concater == "true":
            command = f"""
            source $(conda info --base)/etc/profile.d/conda.sh
            conda activate bactflow
            if [ ! -d {out_dir} ]; then
                mkdir -p {out_dir}
            fi
            {base_dir}/concater.sh -g {fastq_dir} -c {cpus} -e {extension}
            
            """
        else:
            command = f"""
            source $(conda info --base)/etc/profile.d/conda.sh
            conda activate bactflow
            if [ ! -d {out_dir} ]; then
                mkdir -p {out_dir}
            fi
            
            """
        subprocess.run(command, shell=True, text=True, executable="/bin/bash")

       
        if concater == "true":
            command = f"""
            du -sh {fastq_dir}/pooled/*.fastq 
            """
        else:
             command = f"""
            du -sh  {fastq_dir}/*.fastq 
            """

        result = subprocess.run(command, shell=True, text=True, capture_output=True)
        fastq_files = result.stdout.strip().split("\n")
        table_html = """<table class='display table table-striped table-bordered nowrap table-hover' id='read-fastq' border='0.5'>
    <thead>
        <tr>
            <th>Index</th>
            <th>List of Fastq files</th>
            <th>File size</th>
        </tr>
    </thead>
    <tbody>"""
        for i, entry in enumerate(fastq_files, start=1):
            if entry:
                size, filepath = entry.split("\t", 1)
                basename = os.path.basename(filepath)
                table_html += f"<tr><td>{i}</td><td>{basename}</td><td>{size}</td></tr>"
        table_html += "</tbody></table>" 
            
        
        return jsonify({"html_table": table_html, "fastq_files": fastq_files}), 200

# Trim the list
@app.route('/trim-list', methods = ['POST', 'GET'])
def trim_list():
    if request.method == 'POST':
        fastq_dir = request.form.get("fastq_dir")
        concater = request.form.get("concat_reads")
        threshold = request.form.get("size-threshold")

    if concater == "true":
        command = f"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate bactflow
        {base_dir}/read_filter.sh -d {fastq_dir}/pooled -t {threshold} 
        """
    
    else:
        command = f"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate bactflow
        {base_dir}/read_filter.sh -d {fastq_dir} -t {threshold} 
        """
    print(threshold)

    subprocess.run(command, shell=True, text=True, executable="/bin/bash")
    return jsonify({"status": "completed"}), 200


@app.route('/reads-stat', methods = ['POST', 'GET'])
def stat_reads():
    if request.method == 'POST':
        fastq_dir = request.form.get("fastq_dir")
        out_dir = request.form.get("out_dir")
        cpus = request.form.get("cpus")
        concater = request.form.get("concat_reads")
        stats_file = f"{out_dir}/seqkit_stats.tsv"
        
        if concater == "true":
            command = f"""
            source $(conda info --base)/etc/profile.d/conda.sh
            conda activate bactflow
            seqkit stats {fastq_dir}/pooled/*.fastq -a -e -j {cpus} > {out_dir}/seqkit_stats.tsv
            
            """
        else:
            command = f"""
            source $(conda info --base)/etc/profile.d/conda.sh
            conda activate bactflow
            
            seqkit stats {fastq_dir}/*.fastq -a -e  -j {cpus}  > {out_dir}/seqkit_stats.tsv
            
            """
        

        if not os.path.exists(stats_file):
            subprocess.run(command, shell=True, text=True, executable="/bin/bash")
            df = pd.read_csv(f"{out_dir}/fastq_df.csv", sep = "\t")
            if df.empty or df.shape[0] == 0:
                return jsonify({"message": "Stats generation empty", "html_output": ""}), 500
            
        if os.path.exists(stats_file):
            df =  pd.read_csv(stats_file, sep="\s+", engine="python")
            if not df.empty or df.shape[0] != 0:
                html_file = f"{out_dir}/seqkit_stats.html"
                df.to_html(html_file, index=False)

                df['sum_len'] = df['sum_len'].astype(str).str.replace(",", "", regex=True).str.strip().astype(int) 
                df['sum_len'] = (
                df['sum_len'].astype(str)  # Convert to string first
                .str.replace(',', '', regex=True)  # Remove thousands separator
                .str.strip()  # Remove spaces
                .astype(int)  # Convert to integer
            )

           
                df = df.to_dict(orient="records")
            
                styled_html = """ <div id="seqkit_table" class="tabdiv">
                    <div class="table-title">
                        <h3 style="display: flex; margin-bottom: 0;">{{ tabnumber }}</h3>
                        <p>{{ tabcaption }}</p><br>
                    </div>
                    <table id="stats-table" class="display table table-striped table-bordered nowrap table-hover">
                        <thead>
                            <tr>{% for column in table[0].keys() %}<th>{{ column }}</th>{% endfor %}</tr>
                        </thead>
                        <tbody>
                            {% for row in table %}
                            <tr>{% for value in row.values() %}<td>{{ value }}</td>{% endfor %}</tr>
                            {% endfor %}
                        </tbody>
                    </table>
                </div>
                """
                return jsonify({"message": "Stats generated", "html_output": render_template_string(
                    styled_html, id="seqkit_table", tabnumber = "Table 1: Read statistics.", tabcaption="Statistics of raw fastq files were generated by seqkit stats -a command. ", table=df
                )}), 200
            else:
                return jsonify({"message": "Stats table is empty", "html_output": ""}), 200
        else:
            return jsonify({"message": "Stats generation failed", "html_output": ""}), 500
    


# plot quality

def genome_read_quality(fastq_file):
    mean_qualities = []
    length = []
    basename = os.path.basename(fastq_file)
    file_name = os.path.splitext(basename)[0]
    open_func = gzip.open if fastq_file.endswith(".gz") else open 

    with open_func(fastq_file, "rt") as fq:
        for rec in SeqIO.parse(fq, "fastq"):
            #phred = np.mean(rec.letter_annotations["phred_quality"])
            lengths  = len(rec.seq)
            phred_all = np.fromiter(rec.letter_annotations["phred_quality"], dtype=int, count=len(rec.seq))
            mean_err = np.mean(np.power(10, phred_all/-10))
            
            mean_qualities.append((file_name, lengths, -10*np.log10(mean_err)))
            length.append(len(rec.seq))   
    return mean_qualities
# process fastq folder

def process_fq_folder(fastq_folder, threads = 4):
    """process all fastq files in a folder at once"""
    all_quality_data = []
    

    fastq_files = [os.path.join(fastq_folder, f) for f in os.listdir(fastq_folder) if f.endswith((".fastq", ".fq", ".gz", ".fastq.gz", ".fq.gz"))]
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        results = executor.map(genome_read_quality, fastq_files)
    for res in results:
        all_quality_data.extend(res)
    return pd.DataFrame(all_quality_data, columns=["file_name", "read_length", "read_quality"])

@app.route("/plot-qual", methods = ["POST", "GET"])
def plot_qual():
    fastq_dir = request.form.get("fastq_dir")
    out_dir = request.form.get("out_dir")
    cpus = int(request.form.get("cpus", 4))
    concater = request.form.get("concat_reads")
    
    if concater == "true":
        fastq_dir = f"{fastq_dir}/pooled"
    else:
        fastq_dir
    if os.path.exists(os.path.join(out_dir, "fastq_df.csv")):
         
        df = pd.read_csv(f"{out_dir}/fastq_df.csv", sep = "\t")
        if df.empty or df.shape[0] == 0:
            return jsonify({"error": "FASTQ processing resulted in an empty dataset."})
    else:

        df = process_fq_folder(fastq_folder=fastq_dir, threads=cpus)
        if df.empty or df.shape[0] == 0:
            return jsonify({"error": "No data available for visualization."})
        df.to_csv(f"{out_dir}/fastq_df.csv", sep="\t", index=False)

    fig = px.box(data_frame=df, x = "file_name", y = "read_quality", color="file_name", title= "Read quality box plot ")
    fig.update_layout(
        autosize=True,
        height=None,
        width=None,
        # margin=dict(l=5, r=5, t=5, b=5),
        xaxis_title="Fastq files",
        yaxis_title="Read Quality"
    )

    fig_heat = px.density_heatmap(data_frame=df, x = "read_quality", y = "read_length", marginal_x="histogram", marginal_y="histogram", facet_col="file_name", title= "Heatmap density plot of read length vs. quality | Multifaceted")

    fig_heat.update_layout(
        autosize=True,
        height=None,
        width=None,
        # margin=dict(l=5, r=5, t=5, b=5),
        yaxis_title="Read length"
    )

    fig_heat_pool = px.density_heatmap(data_frame=df, x = "read_quality", y = "read_length", title= "Heatmap density plot of read length vs. quality | Pooled", marginal_x="histogram", marginal_y="histogram")

    fig_heat_pool.update_layout(
            autosize=True,
            height=None,
            width=None,
            # margin=dict(l=5, r=5, t=5, b=5),
            xaxis_title="Read quality",
            yaxis_title="Read length"
        )
        

    graphJSON = json.dumps(fig, cls = py.utils.PlotlyJSONEncoder)
    graphJSON2 = json.dumps(fig_heat, cls=py.utils.PlotlyJSONEncoder)
    graphJSON3 = json.dumps(fig_heat_pool, cls=py.utils.PlotlyJSONEncoder)
    return jsonify({"graph": graphJSON, "graph_heat": graphJSON2, "graph_heat2": graphJSON3})

if __name__ == '__main__':
    Timer(1, open_browser).start()

    app.run(debug = True, port = 5000, host = "0.0.0.0", use_reloader = True)

