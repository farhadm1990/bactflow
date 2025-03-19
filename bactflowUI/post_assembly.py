
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
    # 'pycirclize': 'pycirclize',
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


from flask import Flask, render_template, request, redirect, Response, send_from_directory, stream_with_context, jsonify, render_template_string
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
import threading
import base64
from io import BytesIO
import glob


base_dir = os.path.abspath(os.path.dirname(__file__))# we can have access to all files from everywhere
app = Flask(__name__, 
            template_folder = os.path.join(base_dir, "templates"),
            static_folder = os.path.join(base_dir, "static"))





#assembly
@app.route('/',  methods =['GET', 'POST'])
def assembly():

    return render_template('post-assembly.html')


def open_browser():
    webbrowser.open(url="http://127.0.0.1:5001/", new = 2, autoraise=True) # new 2 opens tab while new 1 opens window


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
        text=True
    )
    process_status["pid"]=process.pid
    output_history[:] = []
    for line in iter(process.stdout.readline, ""):
        output_queue.put(line.strip())
        output_history.append(line.strip())

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
            prok_annot = request.form.get('prok_annot', 'false')
            checkm_db = request.form.get('checkm_db', "")
            gtdbtk_data_path = request.form.get('gtdbtk_data_path', "")
            genome_dir = request.form.get("genome_dir")
            run_checkm = request.form.get("run_checkm", "false")
            run_quast = request.form.get("run_quast", "true")
            resume_run = request.form.get("resume_run", "true")
            run_bakta = request.form.get("bakta_annot", "false")
            command = None

            command = f"""if [ ! -d {out_dir} ]; then mkdir -p {out_dir}; fi && 
                    nextflow run {base_dir}/main.nf \
                    --setup_only {setup_only} \
                    --cpus {cpus} \
                    --out_dir {out_dir} \
                    --genome_extension {genome_extension} \
                    --checkm_lineag_check {checkm_lineag_check} \
                    --run_flye false \
                    --circle_genome {circle_genome} \
                    --tax_class {tax_class} \
                    --prok_annot {prok_annot} \
                    --run_checkm {run_checkm} \
                    --checkm_db {checkm_db} \
                    --gtdbtk_data_path {gtdbtk_data_path} \
                    --run_quast {run_quast} \
                    --genome_dir {genome_dir}"""
            if resume_run:
                command = command + " -resume"
            else:
                command = command
                
            output_history[:] = []
            
            back_process = Process(target=run_bact, args=(command, process_status, output_queue, output_history))
            back_process.start()
        
            command = None
            return "Bactflow started successfully!\n", 200
            
        if action == "help":
            command = f"nextflow run {base_dir}/main.nf --help"
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
                print
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
            "prok_annot": True if request.form.get("prok_annot") == "true" else False,
            "run_checkm": True if request.form.get("run_checkm") == "true" else False,
            "run_bakta": True if request.form.get("run_bakta") == "true" else False,
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

# Prokka
@app.route("/check-prok", methods = ["POST"])  
def check_prok():
    out_dir = request.form.get("out_dir")
    gene_type = request.form.get("gene_type", 'CDS')

    if not out_dir:
        return jsonify({"exists":False, "error": "Missing out_dir"}), 400
    prok_path = os.path.join(out_dir, "prokk_out")
    if os.path.exists(prok_path):
        try:

            command = f"""
            {base_dir}/gene_counter_prokka.py -d {prok_path} -t {gene_type} -o {out_dir}/gene_count
            """
            
            gene_count = os.path.join(out_dir, "gene_count.tsv")
            subprocess.run(command, shell=True, text=True, check=True, capture_output=True)
            if not os.path.exists(gene_count):
                subprocess.run(command, shell=True, text=True, check=True, capture_output=True)

            df = pd.read_csv(gene_count, sep="\t")
            
            if os.path.exists(gene_count):
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
                    table_html, id="prokka-tab", tabnumber = "Table 1: Count table of annotated genes by Prokka.", tabcaption="This is the abundance of all genes (selected based on gene type input) at all genomes.", table=table_data
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
        return jsonify({"exists": False, "error": "Prokka annotation directory doesn't exist!"}), 400

# circular
@app.route("/circular", methods=["POST"])
def circular():
    out_dir = request.form.get("out_dir")
    gbk_dir = os.path.join(out_dir, "prokk_out")
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

  
    if os.path.exists(params_file):
        with open(params_file, "r") as f:
            last_params = json.load(f)
    else:
        last_params = {}


    if not os.path.exists(crc_plt) or params != last_params:
        print("🔄 Parameters changed or plot missing. Regenerating plot...")
        
   
        with open(params_file, "w") as f:
            json.dump(params, f, indent=4)

        # source $(conda info --base)/etc/profile.d/conda.sh
        # conda activate bactflow
        command = f"""
        
       
        conda activate bactflow
        {base_dir}/circular_plotter.py -d {gbk_dir} -o {out_dir} \
            {"--add_gc" if params['add_gc'] else ""} \
            {"--add_skew" if params['add_skew'] else ""} \
            --dpi {params["dpi"]} \
            --interval {params["interval"]} \
            --f_color "{params["f_color"]}" \
            --r_color "{params["r_color"]}"
        """
      
        subprocess.run(command, shell=True, check=True)

  
    if os.path.exists(crc_plt):
        with open(crc_plt, "rb") as img_file:
            img_base = base64.b64encode(img_file.read()).decode("utf-8")
        return jsonify({"plot": f"data:image/png;base64,{img_base}"})
    else:
        return jsonify({"error": "Error: Circular plot doesn't exist!"}), 404



    
@app.route("/taxa-report", methods = ["POST"])
def taxa_report():
    out_dir = request.form.get("out_dir")
    tax_file = os.path.join(out_dir, "taxonomy.tsv")

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
    
    output = os.path.join(out_dir, "strain_finder")
    enzyme_file = request.form.get("enzyme_loc")
    count_tab = os.path.join(output, "abundance.tsv")
    plot = os.path.join(output, "requested_genes_abundance.png")

    command = f"""
    {base_dir}/strain_finder.sh -d {gene_files} -g {gene_type} -f abundance -o {output} -e {enzyme_file}  -c false -p {prevalance}

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

    command = f"""
    {base_dir}/strain_finder.sh -d {gene_files} -g {gene_type} -f prevalance -o {output} -e {enzyme_file}  -c false -p {prevalance}

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
