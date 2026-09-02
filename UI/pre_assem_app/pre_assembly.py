#!/usr/bin/env python3

import sys, subprocess, importlib, os, shlex, re, shutil
from flask import Flask, render_template, request, redirect, Response, stream_with_context, jsonify, render_template_string
import plotly.express as px
import concurrent.futures
import plotly as py
from Bio import SeqIO
import json
import gzip
from datetime import datetime, timezone
import pandas as pd
import numpy as np
import os
import subprocess
import sys
from multiprocessing import Process, Manager, Queue
import time
import signal
from threading import Timer
import webbrowser


base_dir = os.path.abspath(os.path.dirname(__file__))
app = Flask(__name__, 
            template_folder = os.path.join(base_dir, "templates"),
            static_folder = os.path.join(base_dir, 'static'))

BACTFLOW_RUNTIME_SH = os.path.join(base_dir, "bactflow_runtime.sh")

NF_JAVA_SETUP = f"""
source "{BACTFLOW_RUNTIME_SH}"
bactflow_prepare_nextflow || exit 1
echo "Using Java for Nextflow: $JAVA_CMD"
"$JAVA_CMD" -version
echo "Using Nextflow: $BACTFLOW_NEXTFLOW_BIN"
"$BACTFLOW_NEXTFLOW_BIN" -version
export NXF_ANSI_LOG=false
"""


def with_nextflow_java(command):
    return NF_JAVA_SETUP + "\n" + command


FASTQ_EXTS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def list_fastq_files(directory):
    if not directory or not os.path.isdir(directory):
        return []
    files = []
    for name in sorted(os.listdir(directory)):
        path = os.path.join(directory, name)
        if os.path.isfile(path) and name.lower().endswith(FASTQ_EXTS):
            files.append(path)
    return files


def stats_error(message, status=400):
    return jsonify({"error": message, "message": message, "html_output": ""}), status


def is_illumina_form(form):
    return (form.get("read_platform") or "").strip().lower() == "illumina"


def concat_enabled(form):
    if is_illumina_form(form):
        return False
    return (form.get("concat_reads") or "false").lower() == "true"


def illumina_filter_enabled(form):
    return is_illumina_form(form) and (form.get("illumina_filter") or "false").lower() == "true"


IN_DOCKER = os.environ.get("BACTFLOW_IN_DOCKER") == "1" or os.path.exists("/.dockerenv")


def conda_shell(command):
    return f"""
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bactflow
{command}
"""


def _strip_ansi(text):
    return re.sub(r"\x1b\[[0-9;]*m", "", text or "")


def run_seqkit(command):
    result = subprocess.run(
        conda_shell(command),
        shell=True,
        text=True,
        executable="/bin/bash",
        capture_output=True,
    )
    if result.returncode == 0:
        return None
    detail = _strip_ansi((result.stderr or result.stdout or "seqkit failed").strip())
    return detail or "seqkit failed"


def pair_illumina_files(fastq_files):
    grouped = {}
    singles = []
    for path in fastq_files:
        mate = illumina_mate_label(path)
        sample = illumina_sample_id(path)
        if mate in ("R1", "R2"):
            grouped.setdefault(sample, {})[mate] = path
        else:
            singles.append(path)
    pairs = []
    for sample, mates in grouped.items():
        if "R1" in mates and "R2" in mates:
            pairs.append((sample, mates["R1"], mates["R2"]))
        else:
            singles.extend(mates.values())
    return pairs, singles


def run_illumina_filter(form):
    fastq_dir = (form.get("fastq_dir") or "").strip()
    out_dir = (form.get("out_dir") or "./bactflow_out").strip()
    try:
        cpus = max(1, int(form.get("cpus") or 1))
    except ValueError:
        cpus = 1
    try:
        min_len = max(1, int(float(form.get("illumina_min_length") or 50)))
    except ValueError:
        min_len = 50
    try:
        min_qual = float(form.get("illumina_min_quality") or 20)
    except ValueError:
        min_qual = 20.0

    if not fastq_dir or not os.path.isdir(fastq_dir):
        return None, "Choose an existing FASTQ directory first."

    os.makedirs(out_dir, exist_ok=True)
    out_dir = os.path.abspath(out_dir)
    dest = os.path.join(out_dir, "illumina_filt")
    os.makedirs(dest, exist_ok=True)
    fastq_files = list_fastq_files(fastq_dir)
    if not fastq_files:
        return None, f"No FASTQ files found in {fastq_dir}."

    fingerprint = f"{_file_fingerprint(fastq_files)}|len={min_len}|q={min_qual}"
    marker = os.path.join(dest, ".filter_meta.json")
    if os.path.isfile(marker):
        try:
            with open(marker, "r", encoding="utf-8") as handle:
                meta = json.load(handle)
            if meta.get("fingerprint") == fingerprint and list_fastq_files(dest):
                return dest, None
        except (OSError, json.JSONDecodeError, ValueError):
            pass

    # seqkit pair -f -O dest deletes dest, so work files must live outside it.
    work = os.path.join(out_dir, ".illumina_filt_work")
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work, exist_ok=True)
    shutil.rmtree(os.path.join(dest, "tmp"), ignore_errors=True)
    for stale in list_fastq_files(dest):
        try:
            os.remove(stale)
        except OSError:
            pass

    pairs, singles = pair_illumina_files(fastq_files)

    def seqkit_seq(src, dst):
        return (
            f"seqkit seq -g --quiet -m {min_len} -Q {min_qual} -j {cpus} "
            f"{shlex.quote(src)} -o {shlex.quote(dst)}"
        )

    def fail(detail):
        shutil.rmtree(work, ignore_errors=True)
        return None, f"Illumina filtering failed: {detail}"

    for sample, r1, r2 in pairs:
        r1_out = os.path.join(work, f"{sample}_R1.fastq.gz")
        r2_out = os.path.join(work, f"{sample}_R2.fastq.gz")
        err = run_seqkit(seqkit_seq(r1, r1_out))
        if err:
            return fail(err)
        err = run_seqkit(seqkit_seq(r2, r2_out))
        if err:
            return fail(err)
        if not os.path.isfile(r1_out) or not os.path.isfile(r2_out):
            return fail(
                f"{sample}: no reads passed the length/quality cutoffs. "
                "Lower min length or min quality."
            )
        pair_dir = os.path.join(work, f"paired_{sample}")
        os.makedirs(pair_dir, exist_ok=True)
        err = run_seqkit(
            "seqkit pair --quiet "
            f"-1 {shlex.quote(r1_out)} -2 {shlex.quote(r2_out)} "
            f"-O {shlex.quote(pair_dir)} -j {cpus}"
        )
        if err:
            return fail(err)
        moved = False
        for name in os.listdir(pair_dir):
            src = os.path.join(pair_dir, name)
            if os.path.isfile(src) and name.lower().endswith(FASTQ_EXTS):
                shutil.move(src, os.path.join(dest, name))
                moved = True
        if not moved:
            return fail(f"{sample}: pairing produced no FASTQ files.")

    for path in singles:
        name = os.path.basename(path)
        if not name.endswith(".gz"):
            name = name + ".gz"
        dst = os.path.join(dest, name)
        err = run_seqkit(seqkit_seq(path, dst))
        if err:
            return fail(err)

    shutil.rmtree(work, ignore_errors=True)

    kept = list_fastq_files(dest)
    if not kept:
        return None, "Illumina filtering removed every read. Lower the length/quality cutoffs."

    with open(marker, "w", encoding="utf-8") as handle:
        json.dump({"fingerprint": fingerprint, "files": [os.path.basename(p) for p in kept]}, handle)
    return dest, None


def prepare_reads_dir(form):
    fastq_dir = (form.get("fastq_dir") or "").strip()
    if not fastq_dir or not os.path.isdir(fastq_dir):
        return None, "Choose an existing FASTQ directory first."
    if illumina_filter_enabled(form):
        return run_illumina_filter(form)
    reads_dir = os.path.join(fastq_dir, "pooled") if concat_enabled(form) else fastq_dir
    if not os.path.isdir(reads_dir):
        return None, (
            "No pooled/ folder yet. Click “Initiate read cat” first, "
            "or set Concat Reads to False."
        )
    return reads_dir, None

manager = Manager()
process_status = manager.dict({"running": False})
output_queue = Queue()

def run_bact(command, process_status, output_queue):
    """To run bactflow in mutliprocess"""
    process_status["running"] = True
    process = subprocess.Popen(
        command,
        shell=True,
        executable="/bin/bash",
        stdout=subprocess.PIPE,
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
    if IN_DOCKER or shutil.which("seqkit"):
        return jsonify({"status": "success", "installed": True}), 200
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

        if error and not output:
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
            


        command = with_nextflow_java(f"""
        
        nextflow run {base_dir}/main.nf --setup_only true --out_dir {out_dir}
        """)
    
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
        concater = "true" if concat_enabled(request.form) else "false"
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
            du -sh {fastq_dir}/pooled/*.fastq* 
            """
        else:
             command = f"""
            du -sh  {fastq_dir}/*.fastq* 
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
        concater = "true" if concat_enabled(request.form) else "false"
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
 

    subprocess.run(command, shell=True, text=True, executable="/bin/bash")
    return jsonify({"status": "completed"}), 200


@app.route("/filter-illumina", methods=["POST"])
def filter_illumina():
    if not is_illumina_form(request.form):
        return jsonify({"error": "Illumina filtering is only available when Read type is Illumina."}), 400
    dest, err = run_illumina_filter(request.form)
    if err:
        return jsonify({"error": err, "message": err}), 400
    kept = list_fastq_files(dest)
    return jsonify({
        "status": "completed",
        "message": f"Filtered Illumina reads written to {dest} ({len(kept)} files).",
        "reads_dir": dest,
    }), 200


@app.route('/reads-stat', methods = ['POST', 'GET'])
def stat_reads():
    if request.method != 'POST':
        return stats_error("POST required", 405)

    try:
        try:
            cpus = max(1, int(request.form.get("cpus") or 1))
        except ValueError:
            cpus = 1

        reads_dir, err = prepare_reads_dir(request.form)
        if err:
            return stats_error(err)

        out_dir = os.path.abspath(request.form.get("out_dir") or "./bactflow_out")
        os.makedirs(out_dir, exist_ok=True)
        stats_file = os.path.join(out_dir, "seqkit_stats.tsv")

        fastq_files = list_fastq_files(reads_dir)
        if not fastq_files:
            return stats_error(
                f"No FASTQ files found in {reads_dir}. "
                "Expected .fastq, .fastq.gz, .fq, or .fq.gz."
            )

        quoted_files = " ".join(shlex.quote(path) for path in fastq_files)
        command = f"""
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate bactflow
        seqkit stats {quoted_files} -a -T -e -j {cpus}
        """
        result = subprocess.run(
            command,
            shell=True,
            text=True,
            executable="/bin/bash",
            capture_output=True,
        )
        if result.returncode != 0 or not (result.stdout or "").strip():
            detail = (result.stderr or result.stdout or "seqkit produced no output").strip()
            if os.path.isfile(stats_file) and os.path.getsize(stats_file) == 0:
                try:
                    os.remove(stats_file)
                except OSError:
                    pass
            return stats_error(f"seqkit stats failed: {detail}", 500)

        with open(stats_file, "w", encoding="utf-8") as handle:
            handle.write(result.stdout)

        try:
            df = pd.read_csv(stats_file, sep="\t")
        except pd.errors.EmptyDataError:
            return stats_error(
                "seqkit wrote an empty table. Check that the FASTQ paths are readable.",
                500,
            )

        if df is None or df.empty:
            return stats_error("seqkit stats table is empty.", 500)

        if "sum_len" in df.columns:
            df["sum_len"] = (
                df["sum_len"].astype(str)
                .str.replace(",", "", regex=False)
                .str.strip()
                .astype(int)
            )

        html_file = os.path.join(out_dir, "seqkit_stats.html")
        df.to_html(html_file, index=False)
        table = df.to_dict(orient="records")
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
        return jsonify({
            "message": "Stats generated",
            "html_output": render_template_string(
                styled_html,
                id="seqkit_table",
                tabnumber="Table 1: Read statistics.",
                tabcaption="Statistics of raw FASTQ files were generated by seqkit stats -a -T.",
                table=table,
            ),
        }), 200
    except Exception as exc:
        return stats_error(str(exc), 500)
 


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
    

    fastq_files = list_fastq_files(fastq_folder)
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        results = executor.map(genome_read_quality, fastq_files)
    for res in results:
        all_quality_data.extend(res)
    return pd.DataFrame(all_quality_data, columns=["file_name", "read_length", "read_quality"])


ILLUMINA_NAME_RE = re.compile(
    r"(?:_R[12]|_r[12]|_[12])(?:_001)?\.(?:fastq|fq)(?:\.gz)?$|illumina",
    re.IGNORECASE,
)
ILLUMINA_MAX_READS = 120000


def _open_fastq(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def illumina_mate_label(filename):
    name = os.path.basename(filename)
    if re.search(r"(_R2|_r2|_2)(?:_001)?\.(?:fastq|fq)", name, re.I):
        return "R2"
    if re.search(r"(_R1|_r1|_1)(?:_001)?\.(?:fastq|fq)", name, re.I):
        return "R1"
    return "SE"


def illumina_sample_id(filename):
    name = os.path.basename(filename)
    name = re.sub(r"\.(fastq|fq)(\.gz)?$", "", name, flags=re.I)
    name = re.sub(r"(_R[12]|_r[12]|_[12])$", "", name)
    return name


def peek_median_length(fastq_file, n=80):
    lengths = []
    try:
        with _open_fastq(fastq_file) as handle:
            for rec in SeqIO.parse(handle, "fastq"):
                lengths.append(len(rec.seq))
                if len(lengths) >= n:
                    break
    except OSError:
        return None
    return float(np.median(lengths)) if lengths else None


def detect_read_platform(fastq_files, declared="auto"):
    declared = (declared or "auto").strip().lower()
    if declared in ("ont", "illumina"):
        return declared
    names = [os.path.basename(path) for path in fastq_files]
    if any(ILLUMINA_NAME_RE.search(name) for name in names):
        return "illumina"
    median = peek_median_length(fastq_files[0])
    if median is not None and median <= 500:
        return "illumina"
    return "ont"


def _file_fingerprint(fastq_files):
    parts = []
    for path in fastq_files:
        try:
            stat = os.stat(path)
            parts.append(f"{os.path.basename(path)}:{stat.st_size}:{int(stat.st_mtime)}")
        except OSError:
            parts.append(os.path.basename(path))
    return "|".join(parts)


def illumina_file_qc(fastq_file, max_reads=ILLUMINA_MAX_READS):
    q_sum = []
    q_n = []
    base_counts = []
    q_hist = np.zeros(46, dtype=np.int64)
    n_reads = 0
    base_idx = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}

    with _open_fastq(fastq_file) as handle:
        for rec in SeqIO.parse(handle, "fastq"):
            seq = str(rec.seq).upper()
            quals = rec.letter_annotations.get("phred_quality") or []
            length = min(len(seq), len(quals))
            if length == 0:
                continue
            if length > len(q_sum):
                extra = length - len(q_sum)
                q_sum.extend([0.0] * extra)
                q_n.extend([0] * extra)
                base_counts.extend([[0, 0, 0, 0, 0] for _ in range(extra)])
            q_total = 0.0
            for i in range(length):
                quality = quals[i]
                q_sum[i] += quality
                q_n[i] += 1
                q_total += quality
                base_counts[i][base_idx.get(seq[i], 4)] += 1
            mean_q = int(round(q_total / length))
            if 0 <= mean_q < len(q_hist):
                q_hist[mean_q] += 1
            n_reads += 1
            if n_reads >= max_reads:
                break

    mean_by_pos = [
        (q_sum[i] / q_n[i]) if q_n[i] else None
        for i in range(len(q_n))
    ]
    composition = []
    for counts in base_counts:
        total = float(sum(counts)) or 1.0
        composition.append({
            "A": 100.0 * counts[0] / total,
            "C": 100.0 * counts[1] / total,
            "G": 100.0 * counts[2] / total,
            "T": 100.0 * counts[3] / total,
            "N": 100.0 * counts[4] / total,
        })
    return {
        "file": os.path.basename(fastq_file),
        "sample": illumina_sample_id(fastq_file),
        "mate": illumina_mate_label(fastq_file),
        "n_reads": n_reads,
        "mean_q_by_pos": mean_by_pos,
        "composition": composition,
        "q_hist": q_hist.tolist(),
    }


def build_illumina_figures(summaries):
    cycle_rows = []
    hist_rows = []
    base_rows = []
    sampled = 0
    for rec in summaries:
        sampled += rec["n_reads"]
        label = f"{rec['sample']} {rec['mate']}".strip()
        for pos, quality in enumerate(rec["mean_q_by_pos"], start=1):
            if quality is None:
                continue
            cycle_rows.append({
                "sample": rec["sample"],
                "mate": rec["mate"],
                "file": label,
                "cycle": pos,
                "mean_quality": round(quality, 2),
            })
        for quality, count in enumerate(rec["q_hist"]):
            if count:
                hist_rows.append({
                    "sample": rec["sample"],
                    "mate": rec["mate"],
                    "file": label,
                    "mean_quality": quality,
                    "reads": int(count),
                })
        for pos, bases in enumerate(rec["composition"], start=1):
            for base, pct in bases.items():
                if base == "N":
                    continue
                base_rows.append({
                    "sample": rec["sample"],
                    "mate": rec["mate"],
                    "file": label,
                    "cycle": pos,
                    "base": base,
                    "percent": round(pct, 2),
                })

    cycle_df = pd.DataFrame(cycle_rows)
    hist_df = pd.DataFrame(hist_rows)
    base_df = pd.DataFrame(base_rows)
    if cycle_df.empty or hist_df.empty or base_df.empty:
        raise ValueError("Illumina QC produced no plottable data. Check that the FASTQ files are readable.")
    caption = f"Sampled up to {ILLUMINA_MAX_READS:,} reads per file ({sampled:,} reads total)."

    fig_cycle = px.line(
        cycle_df,
        x="cycle",
        y="mean_quality",
        color="file",
        line_dash="mate",
        title=f"Illumina per-cycle quality | {caption}",
    )
    fig_cycle.add_hrect(y0=28, y1=42, fillcolor="rgba(18,181,168,0.10)", line_width=0)
    fig_cycle.add_hrect(y0=20, y1=28, fillcolor="rgba(238,184,77,0.12)", line_width=0)
    fig_cycle.add_hrect(y0=0, y1=20, fillcolor="rgba(238,108,77,0.10)", line_width=0)
    fig_cycle.add_hline(y=30, line_dash="dot", line_color="#12b5a8")
    fig_cycle.add_hline(y=20, line_dash="dot", line_color="#ee6c4d")
    fig_cycle.update_layout(
        autosize=True,
        xaxis_title="Cycle",
        yaxis_title="Mean Phred quality",
        yaxis_range=[0, 42],
        legend_title_text="FASTQ",
    )

    fig_hist = px.bar(
        hist_df,
        x="mean_quality",
        y="reads",
        color="file",
        barmode="group",
        title=f"Illumina mean quality distribution | {caption}",
    )
    fig_hist.update_layout(
        autosize=True,
        xaxis_title="Mean read quality",
        yaxis_title="Reads",
        legend_title_text="FASTQ",
    )

    fig_base = px.line(
        base_df,
        x="cycle",
        y="percent",
        color="base",
        facet_col="file",
        facet_col_wrap=3,
        title=f"Illumina per-cycle nucleotide content | {caption}",
        color_discrete_map={"A": "#32a852", "C": "#3b82f6", "G": "#111827", "T": "#ef4444"},
    )
    fig_base.update_layout(
        autosize=True,
        xaxis_title="Cycle",
        yaxis_title="Base percent",
        legend_title_text="Base",
    )
    fig_base.update_yaxes(range=[0, 100])
    return fig_cycle, fig_hist, fig_base


def build_ont_figures(df):
    fig = px.box(
        data_frame=df,
        x="file_name",
        y="read_quality",
        color="file_name",
        title="ONT read quality box plot",
    )
    fig.update_layout(
        autosize=True,
        xaxis_title="FASTQ files",
        yaxis_title="Read quality",
    )
    fig_heat = px.density_heatmap(
        data_frame=df,
        x="read_quality",
        y="read_length",
        marginal_x="histogram",
        marginal_y="histogram",
        facet_col="file_name",
        title="ONT read length vs quality | per sample",
    )
    fig_heat.update_layout(autosize=True, yaxis_title="Read length")
    fig_heat_pool = px.density_heatmap(
        data_frame=df,
        x="read_quality",
        y="read_length",
        title="ONT read length vs quality | pooled",
        marginal_x="histogram",
        marginal_y="histogram",
    )
    fig_heat_pool.update_layout(
        autosize=True,
        xaxis_title="Read quality",
        yaxis_title="Read length",
    )
    return fig, fig_heat, fig_heat_pool


@app.route("/plot-qual", methods = ["POST", "GET"])
def plot_qual():
    if request.method != "POST":
        return jsonify({"error": "POST required"}), 405

    declared = (request.form.get("read_platform") or "auto").lower()
    try:
        cpus = max(1, int(request.form.get("cpus") or 4))
    except ValueError:
        cpus = 4

    reads_dir, err = prepare_reads_dir(request.form)
    if err:
        return jsonify({"error": err}), 400

    out_dir = os.path.abspath(request.form.get("out_dir") or "./bactflow_out")
    os.makedirs(out_dir, exist_ok=True)

    fastq_files = list_fastq_files(reads_dir)
    if not fastq_files:
        return jsonify({"error": f"No FASTQ files found in {reads_dir}."}), 400

    if is_illumina_form(request.form):
        declared = "illumina"

    platform = detect_read_platform(fastq_files, declared)

    try:
        if platform == "illumina":
            cache_path = os.path.join(out_dir, "illumina_qc.json")
            fingerprint = _file_fingerprint(fastq_files)
            summaries = None
            if os.path.isfile(cache_path):
                try:
                    with open(cache_path, "r", encoding="utf-8") as handle:
                        cached = json.load(handle)
                    if cached.get("fingerprint") == fingerprint:
                        summaries = cached.get("summaries")
                except (OSError, json.JSONDecodeError, ValueError):
                    summaries = None
            if not summaries:
                with concurrent.futures.ThreadPoolExecutor(max_workers=cpus) as executor:
                    summaries = list(executor.map(illumina_file_qc, fastq_files))
                with open(cache_path, "w", encoding="utf-8") as handle:
                    json.dump({"fingerprint": fingerprint, "summaries": summaries}, handle)
            if not summaries:
                return jsonify({"error": "No Illumina reads could be parsed."}), 400
            fig, fig_heat, fig_heat_pool = build_illumina_figures(summaries)
        else:
            fingerprint = _file_fingerprint(fastq_files)
            ont_cache = os.path.join(out_dir, "ont_fastq_df.csv")
            ont_meta = os.path.join(out_dir, "ont_fastq_df.meta.json")
            df = None
            if os.path.isfile(ont_cache) and os.path.isfile(ont_meta):
                try:
                    with open(ont_meta, "r", encoding="utf-8") as handle:
                        meta = json.load(handle)
                    if meta.get("fingerprint") == fingerprint:
                        df = pd.read_csv(ont_cache, sep="\t")
                except (OSError, json.JSONDecodeError, ValueError, pd.errors.EmptyDataError):
                    df = None
            if df is None or df.empty:
                df = process_fq_folder(fastq_folder=reads_dir, threads=cpus)
                if df.empty:
                    return jsonify({"error": "No data available for visualization."}), 400
                df.to_csv(ont_cache, sep="\t", index=False)
                with open(ont_meta, "w", encoding="utf-8") as handle:
                    json.dump({"fingerprint": fingerprint}, handle)
            fig, fig_heat, fig_heat_pool = build_ont_figures(df)
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500

    return jsonify({
        "platform": platform,
        "graph": json.dumps(fig, cls=py.utils.PlotlyJSONEncoder),
        "graph_heat": json.dumps(fig_heat, cls=py.utils.PlotlyJSONEncoder),
        "graph_heat2": json.dumps(fig_heat_pool, cls=py.utils.PlotlyJSONEncoder),
    })


if __name__ == '__main__':
    skip_browser = IN_DOCKER or os.environ.get("BACTFLOW_NO_BROWSER") == "1"
    if not skip_browser:
        Timer(1, open_browser).start()
    app.run(
        debug=not skip_browser,
        port=5000,
        host="0.0.0.0",
        use_reloader=False,
        threaded=True,
    )

