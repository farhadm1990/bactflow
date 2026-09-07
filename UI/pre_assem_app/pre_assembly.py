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
import threading
from threading import Timer
import webbrowser

try:
    import psutil
except ImportError:
    psutil = None


base_dir = os.path.abspath(os.path.dirname(__file__))
if base_dir not in sys.path:
    sys.path.insert(0, base_dir)
ASSEM_APP_DIR = os.path.abspath(os.path.join(base_dir, "..", "assem_app"))
if ASSEM_APP_DIR not in sys.path:
    sys.path.append(ASSEM_APP_DIR)

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


from pool_reads import (
    FASTQ_EXTS,
    PoolReadsError,
    human_size,
    inspect_layout,
    list_read_files,
    pool_fastq_dir,
    resolve_fastq_dir,
)


def list_fastq_files(directory, recurse=False):
    return list_read_files(directory, recurse=recurse)


def stats_error(message, status=400):
    return jsonify({"error": message, "message": message, "html_output": ""}), status


def declared_platform(form):
    return (form.get("read_platform") or "auto").strip().lower()


def is_illumina_form(form):
    return declared_platform(form) == "illumina"


def is_pacbio_form(form):
    return declared_platform(form) == "pacbio"


def pacbio_read_kind(form):
    kind = (form.get("pacbio_read_kind") or "hifi").strip().lower()
    if kind in ("hifi", "corr", "clr"):
        return kind
    return "hifi"


def classify_pacbio_reads_message(fastq_files, num_records=5000):
    """Classify PacBio FASTQs and return a user-facing HiFi vs subread message."""
    try:
        from pacbio_read_check import classify_path, user_facing_summary
    except ImportError:
        return (
            "PacBio classifier unavailable. If mean quality is low (&lt;15), "
            "these are likely subreads/CLR — use pacbio-raw in Assembly."
        )
    results = []
    for path in fastq_files[:20]:
        try:
            results.append(classify_path(path, num_records=num_records))
        except Exception as exc:
            results.append({
                "path": str(path),
                "read_class": "unknown",
                "message": f"{os.path.basename(path)}: classification failed ({exc})",
                "reason": str(exc),
            })
    if not results:
        return None
    return user_facing_summary(results)


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
    try:
        fastq_dir, _correction = resolve_fastq_dir(form.get("fastq_dir"))
    except PoolReadsError as exc:
        return None, str(exc)
    if illumina_filter_enabled(form):
        return run_illumina_filter(form)
    layout = inspect_layout(fastq_dir)
    pooled = os.path.join(fastq_dir, "pooled")
    want_concat = concat_enabled(form) or layout.get("needs_concat")
    if want_concat:
        pooled_files = list_fastq_files(pooled)
        if pooled_files:
            return pooled, None
        extra = layout.get("note") or (
            "FASTQ files are in sample subfolders (e.g. TL110_fastq/*.subreads.fastq)."
        )
        return None, (
            f"{extra} Click “Initiate read cat” first, "
            "or set Concat Reads to False if files are already pooled."
        )
    files = list_fastq_files(fastq_dir)
    if files:
        return fastq_dir, None
    if layout.get("needs_concat"):
        extra = layout.get("note") or "FASTQ files are inside sample subfolders."
        return None, extra + ' Click “Initiate read cat” to pool each sample first.'
    nested = list_fastq_files(fastq_dir, recurse=True)
    if nested:
        return None, (
            f"Found {len(nested)} FASTQ file(s) in subfolders of {fastq_dir}. "
            "Set Concat Reads to True and click “Initiate read cat”."
        )
    return None, f"No FASTQ files found in {fastq_dir}."

manager = Manager()
process_status = manager.dict({"running": False})
output_queue = Queue()


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
    _repair_stdio_for_fork()
    try:
        proc.start()
    except BrokenPipeError:
        _repair_stdio_for_fork()
        proc.start()


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
    from bactflow_open_browser import open_bactflow_browser
    open_bactflow_browser(port=5000)

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
        _start_worker_process(back_process)
      
        return jsonify({"message": "Installation started", "running": True}), 200
    elif request.method == 'GET':
        return jsonify({"running": process_status["running"]})

# concat reads and read stat route
@app.route('/ls-fastq', methods = ['POST', 'GET'])
def ls_fastq():
    if request.method == 'POST':
        try:
            fastq_dir, correction = resolve_fastq_dir(request.form.get("fastq_dir"))
        except PoolReadsError as exc:
            return jsonify({"error": str(exc)}), 400
        extension = (request.form.get("extension") or "auto").strip() or "auto"
        out_dir = (request.form.get("out_dir") or "").strip()
        try:
            cpus = max(1, int(request.form.get("cpus") or 1))
        except ValueError:
            cpus = 1
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

        layout = inspect_layout(fastq_dir, extension)
        want_concat = concat_enabled(request.form) or layout.get("needs_concat")
        try:
            if want_concat:
                pooled = pool_fastq_dir(fastq_dir, extension=extension, cpus=cpus)
                files = [str(path) for path in pooled]
            else:
                files = list_fastq_files(fastq_dir)
                if not files:
                    files = list_fastq_files(fastq_dir, recurse=True)
        except PoolReadsError as exc:
            return jsonify({"error": str(exc) or "Concatenation failed."}), 400
        except Exception as exc:
            return jsonify({"error": str(exc)}), 500

        if not files:
            extra = layout.get("note") or "No FASTQ files found."
            return jsonify({"error": extra}), 400

        table_html = """<table class='display table table-striped table-bordered nowrap table-hover' id='read-fastq' border='0.5'>
    <thead>
        <tr>
            <th>Index</th>
            <th>List of Fastq files</th>
            <th>File size</th>
        </tr>
    </thead>
    <tbody>"""
        fastq_entries = []
        for i, filepath in enumerate(files, start=1):
            try:
                size = human_size(os.path.getsize(filepath))
            except OSError:
                size = "?"
            basename = os.path.basename(filepath)
            table_html += f"<tr><td>{i}</td><td>{basename}</td><td>{size}</td></tr>"
            fastq_entries.append(f"{size}\t{filepath}")
        table_html += "</tbody></table>"

        return jsonify({
            "html_table": table_html,
            "fastq_files": fastq_entries,
            "note": " ".join(x for x in (correction, layout.get("note") or "") if x),
            "looks_like_subreads": bool(layout.get("looks_like_subreads")),
            "extension": layout.get("extension") or extension,
            "needs_concat": bool(layout.get("needs_concat")),
        }), 200

# Trim the list
@app.route('/trim-list', methods = ['POST', 'GET'])
def trim_list():
    if request.method == 'POST':
        try:
            fastq_dir, _correction = resolve_fastq_dir(request.form.get("fastq_dir"))
        except PoolReadsError as exc:
            return jsonify({"error": str(exc)}), 400
        concater = "true" if concat_enabled(request.form) else "false"
        threshold = request.form.get("size-threshold")
        target = os.path.join(fastq_dir, "pooled") if concater == "true" else fastq_dir
        command = f"""
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate bactflow
        {shlex.quote(base_dir)}/read_filter.sh -d {shlex.quote(target)} -t {shlex.quote(str(threshold or "0"))}
        """
        subprocess.run(command, shell=True, text=True, executable="/bin/bash")
        return jsonify({"status": "completed"}), 200
    return jsonify({"status": "idle"}), 200


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
    basename = os.path.basename(fastq_file)
    file_name = os.path.splitext(basename)[0]
    if file_name.endswith(".fastq") or file_name.endswith(".fq"):
        file_name = os.path.splitext(file_name)[0]
    open_func = gzip.open if fastq_file.endswith(".gz") else open

    with open_func(fastq_file, "rt") as fq:
        for rec in SeqIO.parse(fq, "fastq"):
            lengths = len(rec.seq)
            phred_all = np.fromiter(
                rec.letter_annotations["phred_quality"], dtype=float, count=len(rec.seq)
            )
            if lengths == 0 or phred_all.size == 0:
                continue
            mean_err = float(np.mean(np.power(10.0, phred_all / -10.0)))
            mean_err = max(mean_err, 1e-16)
            mean_qualities.append((file_name, lengths, -10.0 * np.log10(mean_err)))
    return mean_qualities


def process_fq_folder(fastq_folder, threads=4):
    """Process all FASTQ files in a folder (ONT-style mean Phred quality)."""
    all_quality_data = []
    fastq_files = list_fastq_files(fastq_folder)
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        results = executor.map(genome_read_quality, fastq_files)
    for res in results:
        all_quality_data.extend(res)
    return pd.DataFrame(all_quality_data, columns=["file_name", "read_length", "read_quality"])


_PACBIO_RQ_RE = re.compile(r"(?:^|\s|/)rq:([0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?)", re.I)
_PACBIO_NP_RE = re.compile(r"(?:^|\s|/)np:(\d+)", re.I)
_PACBIO_EC_RE = re.compile(r"(?:^|\s|/)ec:([0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?)", re.I)


def _pacbio_header_qv(description):
    """Convert PacBio predicted accuracy (rq in [0,1]) to Phred-like QV."""
    match = _PACBIO_RQ_RE.search(description or "")
    if not match:
        return None
    rq = float(match.group(1))
    rq = min(max(rq, 0.0), 1.0 - 1e-12)
    return float(-10.0 * np.log10(1.0 - rq))


def genome_pacbio_metrics(fastq_file):
    """Per-read PacBio metrics: length, QV (prefer rq tag), optional np/ec."""
    rows = []
    basename = os.path.basename(fastq_file)
    file_name = os.path.splitext(basename)[0]
    if file_name.endswith(".fastq") or file_name.endswith(".fq"):
        file_name = os.path.splitext(file_name)[0]
    open_func = gzip.open if fastq_file.endswith(".gz") else open

    with open_func(fastq_file, "rt") as fq:
        for rec in SeqIO.parse(fq, "fastq"):
            length = len(rec.seq)
            if length == 0:
                continue
            desc = f"{rec.id} {rec.description}"
            qv = _pacbio_header_qv(desc)
            q_source = "rq"
            if qv is None:
                phred_all = np.fromiter(
                    rec.letter_annotations.get("phred_quality") or [],
                    dtype=float,
                    count=length,
                )
                if phred_all.size == 0:
                    continue
                mean_err = float(np.mean(np.power(10.0, phred_all / -10.0)))
                mean_err = max(mean_err, 1e-16)
                qv = float(-10.0 * np.log10(mean_err))
                q_source = "phred"
            np_match = _PACBIO_NP_RE.search(desc)
            ec_match = _PACBIO_EC_RE.search(desc)
            rows.append({
                "file_name": file_name,
                "read_length": length,
                "read_quality": qv,
                "quality_source": q_source,
                "num_passes": int(np_match.group(1)) if np_match else np.nan,
                "effective_coverage": float(ec_match.group(1)) if ec_match else np.nan,
            })
    return rows


def process_pacbio_folder(fastq_folder, threads=4):
    all_rows = []
    fastq_files = list_fastq_files(fastq_folder)
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        results = executor.map(genome_pacbio_metrics, fastq_files)
    for res in results:
        all_rows.extend(res)
    if not all_rows:
        return pd.DataFrame(
            columns=[
                "file_name",
                "read_length",
                "read_quality",
                "quality_source",
                "num_passes",
                "effective_coverage",
            ]
        )
    return pd.DataFrame(all_rows)


def _read_n50(lengths):
    lengths = np.sort(np.asarray(lengths, dtype=float))[::-1]
    if lengths.size == 0:
        return 0.0
    total = lengths.sum()
    cum = 0.0
    for value in lengths:
        cum += value
        if cum >= total * 0.5:
            return float(value)
    return float(lengths[-1])


def build_pacbio_figures(df, pacbio_kind="hifi"):
    """NanoPlot-style PacBio QC (De Coster et al. / Galaxy long-read QC).

    Always shows length AND quality together:
      1) bivariate length vs mean quality (log length) — NanoPlot length×Qscore
      2) read length histogram with N50
      3) cumulative yield by length
    Quality uses rq→QV when present, otherwise mean Phred from FASTQ QUAL.
    """
    kind = (pacbio_kind or "hifi").lower()
    kind_label = {"hifi": "HiFi/CCS", "corr": "corrected", "clr": "CLR"}.get(kind, "PacBio")
    q_src = "phred"
    if "quality_source" in df.columns and len(df):
        mode = df["quality_source"].dropna()
        if len(mode):
            q_src = str(mode.mode().iloc[0])
    q_label = "Predicted QV (rq)" if q_src == "rq" else "Mean Phred quality"
    label = f"PacBio ({kind_label})"

    plot_df = df.copy()
    plot_df["log10_length"] = np.log10(plot_df["read_length"].clip(lower=1))

    # Downsample for scatter readability (NanoPlot --downsample style).
    if len(plot_df) > 25000:
        plot_df = plot_df.sample(n=25000, random_state=1)

    fig_lq = px.density_heatmap(
        data_frame=plot_df,
        x="read_quality",
        y="log10_length",
        nbinsx=60,
        nbinsy=60,
        color_continuous_scale="Viridis",
        title=f"{label} | read length vs quality (NanoPlot-style)",
        labels={
            "read_quality": q_label,
            "log10_length": "log10(read length bp)",
        },
    )
    fig_lq.update_layout(
        autosize=True,
        xaxis_title=q_label,
        yaxis_title="log10(read length bp)",
        coloraxis_colorbar_title="Reads",
    )
    # Tick labels as actual bp on a log axis of log10 values
    tick_vals = [2, 3, 4, 5]
    fig_lq.update_yaxes(
        tickmode="array",
        tickvals=tick_vals,
        ticktext=[f"1e{v}" for v in tick_vals],
    )
    if kind in ("hifi", "corr") or q_src == "rq":
        for q, color, text in ((20, "#f59e0b", "Q20"), (30, "#16a34a", "Q30")):
            fig_lq.add_vline(
                x=q,
                line_dash="dash",
                line_color=color,
                annotation_text=text,
                annotation_position="top",
            )

    # Length histogram + N50 (NanoPlot LengthHistogram)
    n50 = _read_n50(df["read_length"])
    median_len = float(np.median(df["read_length"]))
    mean_q = float(np.mean(df["read_quality"]))
    median_q = float(np.median(df["read_quality"]))
    fig_hist = px.histogram(
        data_frame=df,
        x="read_length",
        color="file_name",
        nbins=80,
        barmode="overlay",
        opacity=0.65,
        title=(
            f"{label} | read length histogram  "
            f"(N50={int(n50):,}  median={int(median_len):,}  "
            f"meanQ={mean_q:.1f}  medianQ={median_q:.1f})"
        ),
        log_x=True,
    )
    fig_hist.update_layout(
        autosize=True,
        xaxis_title="Read length (bp, log scale)",
        yaxis_title="Number of reads",
        legend_title_text="FASTQ",
    )
    fig_hist.add_vline(
        x=n50,
        line_dash="dash",
        line_color="#dc2626",
        annotation_text=f"N50={int(n50)}",
        annotation_position="top",
    )

    # Cumulative yield by length (NanoPlot YieldByLength)
    yield_rows = []
    for file_name, group in df.groupby("file_name"):
        lengths = np.sort(group["read_length"].to_numpy(dtype=float))[::-1]
        if lengths.size == 0:
            continue
        cum_gb = np.cumsum(lengths) / 1e9
        if lengths.size > 5000:
            idx = np.linspace(0, lengths.size - 1, 5000).astype(int)
            lengths_p, cum_p = lengths[idx], cum_gb[idx]
        else:
            lengths_p, cum_p = lengths, cum_gb
        for length, cum in zip(lengths_p, cum_p):
            yield_rows.append({
                "file_name": file_name,
                "read_length": float(length),
                "cumulative_gb": float(cum),
            })
    yield_df = pd.DataFrame(yield_rows)
    fig_yield = px.line(
        data_frame=yield_df,
        x="read_length",
        y="cumulative_gb",
        color="file_name",
        title=f"{label} | cumulative yield by length",
        log_x=True,
    )
    fig_yield.update_layout(
        autosize=True,
        xaxis_title="Minimum read length (bp, log scale)",
        yaxis_title="Cumulative yield (Gbp)",
        legend_title_text="FASTQ",
    )
    fig_yield.add_vline(
        x=n50,
        line_dash="dot",
        line_color="#dc2626",
        annotation_text=f"N50={int(n50)}",
        annotation_position="top",
    )

    note = (
        f"NanoPlot-style PacBio QC using {q_label}. "
        f"Reads={len(df):,}; N50={int(n50):,} bp; "
        f"median length={int(median_len):,} bp; "
        f"median Q={median_q:.1f}."
    )
    if q_src != "rq" and median_q < 15:
        note += (
            " Low mean Phred is expected for CLR/SRA FASTQs without rq tags; "
            "HiFi BAM/FASTQ with rq: usually sits at Q≥20."
        )
    return fig_lq, fig_hist, fig_yield, "nanoplot", note


def build_long_read_figures(df, platform="ont", pacbio_kind="hifi"):
    """Length×quality plots (same figure set for ONT and PacBio)."""
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
        yaxis_title="Read quality (mean Phred)",
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
        xaxis_title="Read quality (mean Phred)",
        yaxis_title="Read length",
    )
    return fig, fig_heat, fig_heat_pool


def build_ont_figures(df):
    return build_long_read_figures(df, platform="ont")


ILLUMINA_NAME_RE = re.compile(
    r"(?:_R[12]|_r[12]|_[12])(?:_001)?\.(?:fastq|fq)(?:\.gz)?$|illumina",
    re.IGNORECASE,
)
PACBIO_NAME_RE = re.compile(
    r"(?:pacbio|pac.?bio|(?:^|[_\-./])hifi(?:[_\-./]|$)|(?:^|[_\-./])ccs(?:[_\-./]|$)|(?:^|[_\-./])clr(?:[_\-./]|$)|pb0|_pb_|\.pb\.)",
    re.IGNORECASE,
)
ONT_NAME_RE = re.compile(
    r"(?:\bont\b|nanopore|minknow|dorado|guppy|promethion|minion|gridion)",
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
    if declared in ("ont", "illumina", "pacbio"):
        return declared
    names = [os.path.basename(path) for path in fastq_files]
    if any(PACBIO_NAME_RE.search(name) for name in names):
        return "pacbio"
    if any(ILLUMINA_NAME_RE.search(name) for name in names):
        return "illumina"
    if any(ONT_NAME_RE.search(name) for name in names):
        return "ont"
    median = peek_median_length(fastq_files[0])
    if median is not None and median <= 500:
        return "illumina"
    # Long reads without a clear vendor tag default to ONT.
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
    elif is_pacbio_form(request.form):
        declared = "pacbio"

    platform = detect_read_platform(fastq_files, declared)
    pb_kind = pacbio_read_kind(request.form) if platform == "pacbio" else "hifi"
    viz_note = None

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
            # ONT and PacBio share the same length×quality plots
            fingerprint = _file_fingerprint(fastq_files)
            cache_name = "pacbio_fastq_df.csv" if platform == "pacbio" else "ont_fastq_df.csv"
            meta_name = "pacbio_fastq_df.meta.json" if platform == "pacbio" else "ont_fastq_df.meta.json"
            lr_cache = os.path.join(out_dir, cache_name)
            lr_meta = os.path.join(out_dir, meta_name)
            df = None
            if os.path.isfile(lr_cache) and os.path.isfile(lr_meta):
                try:
                    with open(lr_meta, "r", encoding="utf-8") as handle:
                        meta = json.load(handle)
                    if meta.get("fingerprint") == fingerprint:
                        df = pd.read_csv(lr_cache, sep="\t")
                except (OSError, json.JSONDecodeError, ValueError, pd.errors.EmptyDataError):
                    df = None
            if df is None or df.empty:
                df = process_fq_folder(fastq_folder=reads_dir, threads=cpus)
                if df.empty:
                    return jsonify({"error": "No data available for visualization."}), 400
                df.to_csv(lr_cache, sep="\t", index=False)
                with open(lr_meta, "w", encoding="utf-8") as handle:
                    json.dump({"fingerprint": fingerprint, "platform": platform}, handle)
            fig, fig_heat, fig_heat_pool = build_long_read_figures(df, platform="ont")
            if platform == "pacbio":
                viz_note = classify_pacbio_reads_message(fastq_files)
    except Exception as exc:
        return jsonify({"error": str(exc)}), 500

    return jsonify({
        "platform": platform,
        "pacbio_read_kind": pb_kind if platform == "pacbio" else None,
        "pacbio_viz_mode": None,
        "note": viz_note,
        "graph": json.dumps(fig, cls=py.utils.PlotlyJSONEncoder),
        "graph_heat": json.dumps(fig_heat, cls=py.utils.PlotlyJSONEncoder),
        "graph_heat2": json.dumps(fig_heat_pool, cls=py.utils.PlotlyJSONEncoder),
    })


_SYS_STATS = {}
_SYS_STATS_LOCK = threading.Lock()
_SYS_STATS_JSON = os.path.join(base_dir, "static", "sys_stats.json")
_CPU_SAMPLE = {"idle": None, "total": None}


def _read_proc_cpu_times():
    try:
        with open("/proc/stat", "r", encoding="utf-8") as handle:
            parts = handle.readline().split()
        vals = [float(x) for x in parts[1:8]]
        idle = vals[3] + (vals[4] if len(vals) > 4 else 0.0)
        total = sum(vals)
        return idle, total
    except (OSError, ValueError, IndexError):
        return None, None


def _cpu_percent_fallback():
    idle, total = _read_proc_cpu_times()
    if idle is None:
        return 0.0
    prev_idle = _CPU_SAMPLE["idle"]
    prev_total = _CPU_SAMPLE["total"]
    _CPU_SAMPLE["idle"] = idle
    _CPU_SAMPLE["total"] = total
    if prev_idle is None or prev_total is None:
        return 0.0
    didle = idle - prev_idle
    dtotal = total - prev_total
    if dtotal <= 0:
        return 0.0
    return max(0.0, min(100.0, (1.0 - didle / dtotal) * 100.0))


def _cgroup_memory_bases():
    """Candidate cgroup dirs for this process (Docker / cgroup v1+v2)."""
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
    """Parse sizes like 32g, 32768m, 1024k, or raw bytes → int bytes."""
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
    """Memory allotment passed by bactflow.sh (--memory / auto cap)."""
    for key in ("BACTFLOW_DOCKER_MEMORY", "BACTFLOW_MEMORY"):
        limit = _parse_memory_size(os.environ.get(key))
        if limit and limit > 0:
            return limit
    return None


def _cgroup_memory_bytes():
    """Return (used, limit) for this container's memory cgroup when capped."""
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
                # Unlimited / nonsense sentinels (cgroup v1 often uses ~2^63)
                if limit <= 0 or limit >= (1 << 62):
                    continue
                # If "limit" is essentially the whole host, treat as uncapped
                if host_total and limit >= host_total * 0.98:
                    continue
                return used, limit
            except (OSError, ValueError):
                continue
    return None, None


def _docker_aware_memory(fallback_used, fallback_total, fallback_percent):
    """Prefer Docker cgroup / BACTFLOW_DOCKER_MEMORY allotment over host RAM."""
    used_b, limit_b = _cgroup_memory_bytes()
    env_limit = _env_docker_memory_limit()
    if limit_b is None and env_limit:
        limit_b = env_limit
        if used_b is None:
            # Approximate usage from host/process view, capped to allotment
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
    ncpu = 1
    cpu = 0.0
    ram_percent = 0.0
    ram_used_gb = 0.0
    ram_total_gb = 0.0
    source = "fallback"

    if psutil is not None:
        try:
            vm = psutil.virtual_memory()
            ncpu = psutil.cpu_count(logical=True) or 1
            cpu = float(psutil.cpu_percent(interval=None))
            ram_percent = float(vm.percent)
            ram_used_gb = round(vm.used / (1024 ** 3), 2)
            ram_total_gb = round(vm.total / (1024 ** 3), 2)
            source = "psutil"
        except Exception:
            cpu = _cpu_percent_fallback()
            source = "proc"
    else:
        cpu = _cpu_percent_fallback()

    if ram_total_gb <= 0 and psutil is None:
        try:
            meminfo = {}
            with open("/proc/meminfo", "r", encoding="utf-8") as handle:
                for line in handle:
                    key, val = line.split(":")
                    meminfo[key] = int(val.strip().split()[0]) * 1024
            total = meminfo.get("MemTotal", 0)
            avail = meminfo.get("MemAvailable", meminfo.get("MemFree", 0))
            used = max(total - avail, 0)
            ram_total_gb = round(total / (1024 ** 3), 2)
            ram_used_gb = round(used / (1024 ** 3), 2)
            ram_percent = round((used / total) * 100.0, 1) if total else 0.0
        except (OSError, ValueError, ZeroDivisionError):
            pass

    ram_used_gb, ram_total_gb, ram_percent, mem_src = _docker_aware_memory(
        ram_used_gb, ram_total_gb, ram_percent
    )
    if mem_src == "docker":
        source = f"{source}+docker"

    return {
        "cpu_percent": round(cpu, 1),
        "cpu_count": ncpu,
        "ram_percent": round(ram_percent, 1),
        "ram_used_gb": ram_used_gb,
        "ram_total_gb": ram_total_gb,
        "ram_scope": mem_src,
        "running": False,
        "job_cores": None,
        "job_cpu_percent": None,
        "job_rss_gb": None,
        "job_procs": 0,
        "source": source,
    }


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
            "error": str(exc),
        }


def _sys_stats_loop():
    if psutil is not None:
        try:
            psutil.cpu_percent(interval=None)
        except Exception:
            pass
    _cpu_percent_fallback()
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


@app.route("/bactflow_status", methods=["GET"])
def bactflow_status():
    payload = current_sys_stats()
    payload["status"] = "ok"
    return jsonify(payload)


if __name__ == '__main__':
    # Docker: open via host browser hook. Local: webbrowser/xdg-open.
    # bactflow.sh sets BACTFLOW_NO_BROWSER=1 and opens the tab itself.
    want_browser = (
        os.environ.get("BACTFLOW_NO_BROWSER") != "1"
        or os.environ.get("BACTFLOW_FORCE_BROWSER") == "1"
    )
    if want_browser:
        Timer(1, open_browser).start()
    flask_debug = os.environ.get("BACTFLOW_FLASK_DEBUG") == "1"
    app.run(
        debug=flask_debug,
        port=5000,
        host="0.0.0.0",
        use_reloader=False,
        threaded=True,
    )

