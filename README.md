<div style="text-align: center; margin-top: 0;">
    <img src="https://github.com/farhadm1990/bactflow/blob/main/logo/logo.png" alt="BactFlow Logo" width="400" height="400"/>
    <p><strong>BactFlow</strong> logo was designed by DALL·E :)</p>
</div>



# BactFlow

## Introduction

BactFlow is a workflow with a user-friendly UI for bacterial genome assembly of single-isolate sequencing reads from **Oxford Nanopore (ONT)**, **Illumina**, **PacBio**, and **ONT/PacBio+Illumina hybrid** data. It is built with Nextflow DSL 2 and reads the usual outputs of Guppy and Dorado basecallers.

One UI on **port 5000** holds every page:

| Page | What it does |
| --- | --- |
| **QC** | Concatenate barcode folders, drop tiny files, seqkit stats, read plots |
| **Assembly** | Flye (ONT / PacBio), SPAdes (Illumina), Unicycler (hybrid), optional Medaka / circulate / QUAST |
| **Post-assembly** | QUAST, Bakta, GTDB-Tk, CheckM, MLST, geNomad, SNP finder, Medaka variant calling, genome plots, Generate report (HTML/PDF) |

Switch pages from the top bar. Only **one assembler** runs per assembly invocation. You can reuse the same `--out_dir` so Illumina, ONT, and PacBio FASTAs pile up in the same folder.

## How you run BactFlow

**Every start goes through `./bactflow.sh run`.** That script picks the image or Conda env, mounts your work directory, waits until the UI is up, and opens the browser. Do **not** run `docker run`, `singularity` / `apptainer`, `conda activate` + `python …/app.py`, or `scripts/setup_bactflow_env.sh` by hand.

### Docker image size (multi-stage)

The Dockerfile is **multi-stage**: conda envs, JRE, BLAST, and UI source are built in separate stages; the final image only copies runtime trees (no micromamba, no JDK, no package caches).

The full unified image is still large mainly because of **geNomad → TensorFlow + MMseqs2**, **Medaka → PyTorch**, and PubMLST schemes. Stripping drops unused Qt/PySide, duplicate MLST BLAST/JVM, TensorBoard, and docs — **without removing app features**.

To go roughly **half size**, omit geNomad only:

```sh
BACTFLOW_PROFILE=slim ./bactflow.sh run "${WORK_DIR}" --service docker --local --rebuild
# or:
docker build -f UI/bactflow_app/Dockerfile \
  --build-arg BACTFLOW_PROFILE=slim \
  -t bactflow:v1.0 .
```

| Profile | geNomad / plasmids | Typical size after rebuild |
| --- | --- | --- |
| `full` (default) | yes | smaller than Hub ~7GB; still multi‑GB (TF + Torch) |
| `slim` | no | ~half of full |

```sh
./bactflow.sh run <work_dir> --service docker|apptainer|conda
```

Leave that terminal running. Stop with **Ctrl+C**. Then open http://127.0.0.1:5000.

You only pick **one** `--service`.

| If this is you | Use | Do you need sudo / a Docker group? |
| --- | --- | --- |
| Linux desktop, Docker already works | **`--service docker`** (this is the default) | Yes, once: your user must be in the `docker` group |
| Company / HPC that will **not** give Docker or sudo | **`--service apptainer`** | No. IT may already have Apptainer/Singularity |
| No containers at all; you can install Conda in your home folder | **`--service conda`** | No |

---

# Step 0 — clone the repo (everyone does this once)

```sh
git clone https://github.com/farhadm1990/bactflow.git
cd bactflow
```

Make a work folder. This is where reads and results should live. Use an **absolute** path (it must start with `/`):

```sh
WORK_DIR="$(pwd)"
mkdir -p "${WORK_DIR}/bactflow_out"
echo "Your work directory is: ${WORK_DIR}"
```

Keep that `WORK_DIR` value. You will paste it into every command below.

Then jump to **exactly one** of: [Docker](#1-docker), [Apptainer](#2-apptainer), or [Conda](#3-conda). Each of those sections only shows `./bactflow.sh`.

---

# 1. Docker

Use this if `docker run hello-world` already works on your account. You still start BactFlow with `./bactflow.sh`, not with `docker run`.

### One-time: Docker group (needs sudo, only once)

```sh
sudo groupadd docker
sudo usermod -aG docker $USER
newgrp docker
docker run hello-world
```

If `hello-world` prints a hello message, Docker is ready. Log out and back in if `newgrp` is not enough.

### Start BactFlow

From the `bactflow` folder:

```sh
cd /path/to/bactflow
WORK_DIR="$(pwd)"

./bactflow.sh run "${WORK_DIR}" --service docker --cpus 10 --memory 16g
```

`--service docker` can be omitted; Docker is the default.

Wait until the terminal prints **UI is ready**. A browser tab should open. If it does not, open http://127.0.0.1:5000.

**Path rule (Docker and Apptainer):** the container only sees files **under** `WORK_DIR`. Put FASTQs, genomes, databases (or symlinks to them), and the output folder inside that directory. In the form, use absolute paths such as `${WORK_DIR}/test_reads/ont` and `${WORK_DIR}/bactflow_out`.

First run downloads images from Docker Hub. That can take several minutes.

---

# 2. Apptainer

Use this on HPC or a locked-down Linux box where you **cannot** join the Docker group. You still start BactFlow with `./bactflow.sh run … --service apptainer`. Do not call `singularity` / `apptainer` yourself.

Apptainer (or Singularity) does **not** need sudo or the Docker group. It uses the **same** Hub images as Docker (`farhadm1990/bactflow_*`). `bactflow.sh` converts each image once under `~/.bactflow/sif/` (a `.sif` file, or a sandbox directory if your Singularity cannot mount squashfs) and reuses it. If you want no containers at all, use [Conda](#3-conda).

You need `apptainer` **or** `singularity` on your PATH (`apptainer --version` or `singularity --version`). If both are missing, ask IT to install Apptainer, or use [Conda](#3-conda) instead.

### Start BactFlow

From the `bactflow` folder:

```sh
cd /path/to/bactflow
WORK_DIR="$(pwd)"

./bactflow.sh run "${WORK_DIR}" --service apptainer --cpus 10 --memory 16g
```

`--service singularity` means the same thing as `--service apptainer`.

The first run converts the Hub image. That is slow once. Later starts reuse `~/.bactflow/sif/`. To refresh:

```sh
./bactflow.sh run "${WORK_DIR}" --service apptainer --pull
```

Then open http://127.0.0.1:5000.

Same **path rule** as Docker: keep reads, databases, and `--out_dir` under `WORK_DIR`.

To store converted images somewhere else (for example a scratch disk):

```sh
export BACTFLOW_SIF_DIR=/scratch/$USER/bactflow-sif
./bactflow.sh run "${WORK_DIR}" --service apptainer
```

---

# 3. Conda

Use this when you do not want containers. Tools install into a Conda env named `bactflow` in your account. **No sudo.** The root `config.yml` has the tools for the whole UI. Create the env and start with `./bactflow.sh run` — do not `conda activate` and run the Flask app yourself.

### One-time: Conda itself (skip if `conda` already works)

Miniforge in your home folder, no root:

```sh
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init bash
exec bash
conda --version
```

### Start BactFlow

From the `bactflow` folder. The **first** launch creates the env (this can take a long time). Later launches reuse it.

```sh
cd /path/to/bactflow
WORK_DIR="$(pwd)"

./bactflow.sh run "${WORK_DIR}" --service conda
```

Then open http://127.0.0.1:5000.

To rebuild/update the env from `config.yml` (Java 17–24 and Nextflow 24.10.2 are included):

```sh
./bactflow.sh run "${WORK_DIR}" --service conda --rebuild
```

Env create **does not download** Bakta, GTDB-Tk, CheckM, or geNomad databases. Install those yourself when you need them and paste the paths in the UI.

Conda runs on the host, so it can see **any** path you can see. Absolute paths are still the least confusing in the forms.

---

# Common flags (all services)

```sh
./bactflow.sh run "${WORK_DIR}" --service docker --cpus 16 --memory 32g
./bactflow.sh run "${WORK_DIR}" --service docker --local --rebuild
./bactflow.sh run "${WORK_DIR}" --service apptainer --tag v1.0
./bactflow.sh run "${WORK_DIR}" --service conda --no-browser
```

| Flag | What it does |
| --- | --- |
| `--service NAME` | `docker` (default), `apptainer`, or `conda` |
| `--cpus N` / `--memory SIZE` | Limits for Docker; Apptainer uses them when the runtime supports it |
| `--port PORT` | Host port if 5000 is already taken |
| `--no-browser` | Do not open a tab; copy the URL from the terminal |
| `--local` | Build from this repo’s Dockerfile (docker / apptainer) |
| `--rebuild` | Rebuild the local Docker image, rebuild the converted Apptainer image, or update the Conda env |
| `--tag TAG` | Pin a Hub tag, for example `v1.0` |
| `--pull` | Re-download the Hub image / rebuild the converted Apptainer image |

If the browser does not open, the terminal still prints the URL. If port 5000 is busy, pass `--port` (for example `--port 5010`).

---

# If something failed

| You saw | Likely fix |
| --- | --- |
| I ran `docker run` / `singularity` / `python …app.py` | Stop that process. Start with `./bactflow.sh run "${WORK_DIR}" --service …` |
| Flask prints several `Running on http://…` links | Use the single **Open this URL** / **Access URL** line. Raw `docker run` without `-p` only reaches the container IP (e.g. `http://172.17.0.2:5000`). Prefer `./bactflow.sh run … --service docker`, which maps the port and uses http://127.0.0.1:5000 |
| `python: can't open file '/app/app.py'` (Apptainer) | Update to the latest `bactflow.sh` (it starts `/app/UI/bactflow_app/app.py`). Then `./bactflow.sh run "${WORK_DIR}" --service apptainer --pull` |
| `Docker is not installed` / `daemon is not reachable` | Use `--service apptainer` or `--service conda`, or fix the Docker group |
| `Apptainer/Singularity is not installed` | Install Apptainer, or use `--service conda` |
| `unable to retrieve auth token: invalid username/password` | Stale Docker Hub login in `~/.docker/config.json`. Retry; BactFlow now prefers `docker save` → `.sif` so Singularity never talks to Hub. |
| `client version 1.22 is too old` | Old Singularity talking to a new Docker daemon. Retry: conversion now uses `docker save` (a tar file), not `docker-daemon://`. |
| `Usage: singularity … run` then the UI never starts | Singularity 3.0 has no `--env`. Retry; BactFlow now uses `SINGULARITYENV_*` and `exec`. |
| `failed to mount squashfs filesystem: invalid argument` | Singularity 3.0 cannot mount a `.sif` on a new kernel. Retry; BactFlow extracts a sandbox directory under `~/.bactflow/sif/`. |
| `conda/mamba was not found` | Install Miniforge in `$HOME` (Conda section above) |
| `Port 5000 is already in use` | `./bactflow.sh run "${WORK_DIR}" --service … --port 5010` |
| Form cannot see my FASTQs (docker / apptainer) | Paths must sit **inside** `WORK_DIR` |
| First conda/apptainer start looks “stuck” | First image/env build is slow; wait for **UI is ready** |

---

# Test reads (ONT, Illumina, PacBio)

The repository includes small **demo FASTQs for three technologies** so you can try the UIs without your own data:

```text
test_reads/
├── ont/                         # ONT (Flye / Unicycler long reads)
│   ├── TL110_Native_R1041.fastq.gz
│   ├── TL19_Native_R1041.fastq.gz
│   └── TL29_Native_R1041.fastq.gz
├── illumina/                    # Illumina paired-end (SPAdes / Unicycler short)
│   ├── TL110_Illumina_R1.fastq.gz
│   ├── TL110_Illumina_R2.fastq.gz
│   ├── TL19_Illumina_R1.fastq.gz
│   ├── TL19_Illumina_R2.fastq.gz
│   ├── TL29_Illumina_R1.fastq.gz
│   └── TL29_Illumina_R2.fastq.gz
└── pacbio/                      # PacBio (Flye; HiFi or set pacbio-raw for subreads)
    ├── TL110.fastq.gz
    ├── TL19.fastq.gz
    └── TL29.fastq.gz
```

Sample IDs (`TL110`, `TL19`, `TL29`) match across technologies so hybrid Unicycler can pair ONT + Illumina.

In the UI, paste **absolute** paths (Docker/Apptainer: they must sit under `WORK_DIR`):

| What you are doing | Path to paste |
| --- | --- |
| Pre-assembly / Assembly ONT | `${WORK_DIR}/test_reads/ont` |
| Assembly Illumina | `${WORK_DIR}/test_reads/illumina` |
| Assembly PacBio | `${WORK_DIR}/test_reads/pacbio` |
| Hybrid Unicycler | long: `${WORK_DIR}/test_reads/ont` · short: `${WORK_DIR}/test_reads/illumina` |
| Output directory | `${WORK_DIR}/bactflow_out` |

---

# What the UIs look like

<div style="text-align: center; margin-top: 10;">
    <img src="https://github.com/farhadm1990/bactflow/blob/main/pix/preassem.png" alt="preassem" style="max-width: 100%; height: auto;"/>
    <p><strong>Fig 1. </strong> Pre-assembly module. </p>
</div>
<br>

<div style="text-align: center; margin-top: 10;">
    <img src="https://github.com/farhadm1990/bactflow/blob/main/pix/assem.png" alt="Assem" style="max-width: 100%; height: auto;"/>
    <p><strong>Fig 2. </strong> Assembly module. </p>
</div>
<br>

<div style="text-align: center; margin-top: 10;">
    <img src="https://github.com/farhadm1990/bactflow/blob/main/pix/postassem.png" alt="Post assem" style="max-width: 100%; height: auto;"/>
    <p><strong>Fig 3. </strong> Post-assembly module. </p>
</div>

The assembly image / conda env includes Flye, SPAdes, Unicycler, Circlator, and QUAST. Large taxonomy databases (Bakta, GTDB-Tk, CheckM) are pointed at from the UI; they are not inside the image.

---

# Input layouts

Use **absolute paths**. One file (or one Illumina pair) per sample.

### ONT (Flye or Unicycler long reads)

Already one FASTQ per sample (`--concat_reads false`) — matches `test_reads/ont`:

```text
ont_reads/
├── TL110_Native_R1041.fastq.gz
├── TL19_Native_R1041.fastq.gz
└── TL29_Native_R1041.fastq.gz
```

Guppy/Dorado barcode folders, not yet pooled (`--concat_reads true`). Each subdirectory is one sample; chunks are concatenated into `ont_reads/pooled/`:

```text
ont_reads/
├── barcode01/*.fastq.gz
├── barcode02/*.fastq.gz
└── barcode03/*.fastq.gz
```

### Illumina paired-end (SPAdes, or Unicycler short reads)

R1/R2 names must match. Accepted patterns include `*_R1.fastq.gz` / `*_R2.fastq.gz`, `*_1.fastq.gz` / `*_2.fastq.gz`, and `sample_Illumina_R1.fastq.gz`.

```text
illumina/
├── TL110_Illumina_R1.fastq.gz
├── TL110_Illumina_R2.fastq.gz
├── TL19_Illumina_R1.fastq.gz
└── TL19_Illumina_R2.fastq.gz
```

For hybrid Unicycler, the ONT filename is peeled down to a sample id (`TL110`) and matched to Illumina files such as `TL110_*_R1.fastq.gz`.

### PacBio (Flye)

```text
pacbio/
├── TL110.fastq.gz
├── TL19.fastq.gz
└── TL29.fastq.gz
```

Assembly detects **HiFi vs subread**-style PacBio inputs. Use `--pacbio_read_type pacbio-hifi` (default) for HiFi; if the check reports **subreads**, set **`pacbio-raw`** in the UI / CLI (BactFlow does not convert subreads→HiFi with CCS).

---

# Running from the terminal

This is an optional **no-UI** Nextflow path for people who already have the Conda env from `./bactflow.sh run … --service conda`. To open the UI, use `./bactflow.sh run` only (see [How you run BactFlow](#how-you-run-bactflow)).

```sh
conda activate bactflow
cd /path/to/bactflow
```

Help for the current assembly workflow:

```sh
nextflow run UI/assem_app/main.nf --help
```

Add `-resume` to continue an interrupted run. Use `--cpus` for threads. Set `--out_dir` to an absolute path if you can.

In the **Assembly** and **Post-assembly** UI, **Clean up work dir** (True/False) deletes `out_dir/.nextflow-work` and other temp files after a **successful** run to free disk. Leave it **False** if you plan to use **Resume** on the next run (Resume needs that work directory). Failed or stopped runs keep the work dir so you can debug or resume.

## 1. Pre-assembly

The pre-assembly UI is the usual way to concatenate, trim tiny files, and plot reads. From the terminal you can concatenate barcode folders the same way the UI does:

```sh
bash UI/pre_assem_app/concater.sh \
  -g /home/user/work_dir/ont_reads \
  -c 8 \
  -e .fastq.gz
```

That writes one FASTQ per barcode into `/home/user/work_dir/ont_reads/pooled/`. Then generate seqkit stats:

```sh
mkdir -p bactflow_out
seqkit stats /home/user/work_dir/ont_reads/pooled/*.fastq -a -e -j 8 > bactflow_out/seqkit_stats.tsv
```

## 2. Assembly — examples by read type

Run **one** assembler per command. Reuse `bactflow_out` if you assemble the same samples with more than one tool; FASTA names are tagged (`_flye`, `_spades`, `_unicycler`, `_pacbio`) so they are not overwritten.

Paths below use the repo’s `test_reads/` — replace with your own absolute paths as needed.

### ONT — Flye

Raw ONT (default `--ont_read_type nano-raw`). Use `nano-hq` for Q20+/Dorado HAC/SUP, or `nano-corr` for already-corrected reads.

```sh
nextflow run UI/assem_app/main.nf \
  --run_flye true \
  --run_spades false \
  --run_unicycler false \
  --run_pacbio false \
  --fastq_dir "$(pwd)/test_reads/ont" \
  --concat_reads false \
  --ont_read_type nano-raw \
  --nanofilter true \
  --min_length 1000 \
  --min_quality 16 \
  --circle_genome true \
  --run_quast true \
  --cpus 10 \
  --out_dir "$(pwd)/bactflow_out" \
  -resume
```

Optional Medaka polish (needs the basecaller model that matches the reads):

```sh
nextflow run UI/assem_app/main.nf \
  --run_flye true \
  --fastq_dir "$(pwd)/test_reads/ont" \
  --concat_reads false \
  --ont_read_type nano-hq \
  --medaka_polish true \
  --basecaller_model r1041_e82_400bps_hac_v4.2.0 \
  --circle_genome true \
  --run_quast true \
  --cpus 10 \
  --out_dir "$(pwd)/bactflow_out"
```

Coverage downsample before Flye:

```sh
nextflow run UI/assem_app/main.nf \
  --run_flye true \
  --fastq_dir "$(pwd)/test_reads/ont" \
  --concat_reads false \
  --coverage_filter true \
  --coverage 50 \
  --genome_size 6 \
  --out_dir "$(pwd)/bactflow_out" \
  --cpus 10
```

### Illumina — SPAdes isolate

`--fastq_dir` must be the Illumina paired-end folder (not ONT).

```sh
nextflow run UI/assem_app/main.nf \
  --run_flye false \
  --run_spades true \
  --run_unicycler false \
  --run_pacbio false \
  --fastq_dir "$(pwd)/test_reads/illumina" \
  --circle_genome true \
  --run_quast true \
  --cpus 10 \
  --out_dir "$(pwd)/bactflow_out" \
  -resume
```

### Hybrid ONT + Illumina — Unicycler

Long reads in `--fastq_dir`, short-read pairs in `--short_read_dir`. Sample prefixes must match (for example ONT `TL110_Native_R1041.fastq.gz` with Illumina `TL110_Illumina_R1.fastq.gz`).

```sh
nextflow run UI/assem_app/main.nf \
  --run_flye false \
  --run_spades false \
  --run_unicycler true \
  --run_pacbio false \
  --fastq_dir "$(pwd)/test_reads/ont" \
  --short_read_dir "$(pwd)/test_reads/illumina" \
  --concat_reads false \
  --circle_genome true \
  --run_quast true \
  --cpus 10 \
  --out_dir "$(pwd)/bactflow_out" \
  -resume
```

### PacBio — Flye

`--pacbio_read_type` is `pacbio-hifi` (default), `pacbio-raw`, or `pacbio-corr`.

```sh
nextflow run UI/assem_app/main.nf \
  --run_flye false \
  --run_spades false \
  --run_unicycler false \
  --run_pacbio true \
  --fastq_dir "$(pwd)/test_reads/pacbio" \
  --concat_reads false \
  --pacbio_read_type pacbio-hifi \
  --circle_genome true \
  --run_quast true \
  --cpus 10 \
  --out_dir "$(pwd)/bactflow_out" \
  -resume
```

### Assembly plus taxonomy / annotation / CheckM

Point databases at absolute paths (or at symlinks under the Docker work directory).

```sh
nextflow run UI/assem_app/main.nf \
  --run_flye true \
  --fastq_dir "$(pwd)/test_reads/ont" \
  --concat_reads false \
  --circle_genome true \
  --run_quast true \
  --bakta_annot true \
  --bakta_db /home/user/work_dir/bakta_db \
  --tax_class true \
  --gtdbtk_data_path /home/user/work_dir/gtdbtk_db/release220 \
  --run_checkm true \
  --checkm_db /home/user/work_dir/checkm_db \
  --cpus 10 \
  --out_dir "$(pwd)/bactflow_out"
```

## 3. Post-assembly (already have FASTA files)

Skip assemblers and score genomes in `--genome_dir` (usually `bactflow_out/asm_out_dir/fastas` or `bactflow_out/circulated_fasta`).

```sh
nextflow run UI/post_assem_app/main.nf \
  --run_flye false \
  --genome_dir "$(pwd)/bactflow_out/asm_out_dir/fastas" \
  --genome_extension fasta \
  --out_dir "$(pwd)/bactflow_out" \
  --run_quast true \
  --bakta_annot true \
  --bakta_db /home/user/work_dir/bakta_db \
  --tax_class true \
  --gtdbtk_data_path /home/user/work_dir/gtdbtk_db/release220 \
  --run_checkm true \
  --checkm_db /home/user/work_dir/checkm_db \
  --cpus 10 \
  -resume
```

### Post-assembly UI extras

From **Post-assembly** in the UI you can also:

- **SNP finder** — BWA + bcftools; interactive VCF tables per genome under `out_dir/snps/`
- **Variant calling (Medaka)** — tables from each genome’s **`medaka.vcf`** under `out_dir/vcs/`
- **Circular plot** — from Bakta GBK/GFF outputs
- **Strain finder** — abundance / prevalence tables and plots from enzyme + annotation inputs
- **Clean up work dir** — same as Assembly: remove `.nextflow-work` / temp files after success (keep False for Resume)

---

# Output directory (`bactflow_out`)

Published results land here. Nextflow work files go to `bactflow_out/.nextflow-work`. Set **Clean up work dir = True** in the Assembly or Post-assembly UI to delete that folder (and other temp files) after a successful run.

Folders that only appear when the matching flag is on are marked *(optional)*.

```text
bactflow_out/
├── environment_created
├── seqkit_stats.tsv                          # pre-assembly seqkit table (optional)
├── seqkit_stats.html
│
├── asm_out_dir/
│   └── fastas/                               # all assembler FASTAs, names are kept
│       ├── TL110_flye.fasta
│       ├── TL110_spades.fasta
│       ├── TL110_unicycler.fasta
│       ├── TL110_pacbio.fasta
│       ├── TL19_flye.fasta
│       └── TL19_spades.fasta
│
├── circulated_fasta/                         # --circle_genome true (optional)
│   ├── TL110_flye.fasta
│   ├── TL110_spades.fasta
│   └── TL19_flye.fasta
│
├── quast_stat/                               # --run_quast true (optional)
│   ├── report.html                           # open this in a browser
│   ├── report.tsv
│   ├── report.txt
│   ├── report.pdf
│   ├── icarus.html
│   ├── icarus_viewers/
│   │   └── contig_size_viewer.html
│   └── basic_stats/
│       ├── Nx_plot.pdf
│       ├── GC_content_plot.pdf
│       └── cumulative_plot.pdf
│
├── bakta_out/                                # --bakta_annot true (optional)
│   └── TL110_flye_bakta/
│       ├── TL110_flye.gbff
│       ├── TL110_flye.gff3
│       ├── TL110_flye.faa
│       └── TL110_flye.ffn
│
├── gtdbtk_out/                               # --tax_class true (optional)
│   └── ...                                   # GTDB-Tk classify_wf results
│
├── checkm_out/                               # --run_checkm true (optional)
│   ├── checkm_lineage.txt
│   ├── taxon_tree.newick
│   ├── genome_tree.newick
│   └── genome_tree.tree
│
├── snps/                                     # post-assembly SNP finder (optional)
│   └── <genome>/… .vcf
│
├── vcs/                                      # post-assembly Medaka VC (optional)
│   └── <genome>_vs_<ref>/medaka.vcf
│
└── strain_finder/                            # post-assembly strain finder (optional)
    ├── abundance.tsv
    ├── prevalance.tsv
    ├── requested_genes_abundance.jpeg
    └── requested_genes_prevalence.jpeg
```

What to pick up first:

- **Assemblies:** `bactflow_out/asm_out_dir/fastas/`
- **Start-fixed (circulated) genomes:** `bactflow_out/circulated_fasta/`
- **Assembly QC:** `bactflow_out/quast_stat/report.html`

QUAST is rebuilt from **every** FASTA already in `asm_out_dir/fastas` (or `circulated_fasta` when circulate is on). So SPAdes then Flye into the same `bactflow_out` produces one report with both.

---

# Databases

### GTDB-Tk

Set `--tax_class true` and pass `--gtdbtk_data_path` to the extracted release directory (for example `release220`).

```sh
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz
tar xvzf gtdbtk_data.tar.gz
```

You can also put a symlink to that database inside the Docker work directory.

### CheckM

Set `--run_checkm true` and `--checkm_db` to the CheckM data root.

### Bakta

Set `--bakta_annot true` and `--bakta_db` to a Bakta database directory.

---

# Generic outputs

## Report maker

<div style="text-align: center; margin-top: 10;">
    <img src="https://github.com/farhadm1990/bactflow/blob/main/pix/pan_genome.png" alt="Pangenome" style="max-width: 100%; height: auto;"/>
    <p><strong>Fig 4. </strong> Pangenome of assembled genomes extracted from BactFlow. </p>
</div><br>

<div style="text-align: center; margin-top: 10;">
    <img src="https://github.com/farhadm1990/bactflow/blob/main/pix/tree_genome.png" alt="Phylogenetic tree" style="max-width: 100%; height: auto;"/>
    <p><strong>Fig 5. </strong> Phylogenetic tree of assembled genomes based on their ANI. </p>
</div>

---

## Potential issues

### Browser does not open (Docker)

- Start with `./bactflow.sh run …` (opens the host browser after the UI is ready).
- If the tab still does not appear, copy the URL printed in the terminal.
- Rebuild the local image after UI changes: `./bactflow.sh run <work_dir> --local --rebuild`.

### Paths outside the Docker work directory

Anything not under the absolute `work_dir` passed to `bactflow.sh` is invisible inside the container. Symlink large databases into that tree if needed.

### Conda unbound variable on first Nextflow env setup

On the first launch, process `envSetUP` creates a conda environment called `bactflow`. In some conda setups this error can appear:

```sh
miniconda3/envs/bactflow/etc/conda/deactivate.d/libxml2_deactivate.sh: line 3: xml_catalog_files_libxml2: unbound variable
```

Edit that file as follows and rerun BactFlow:

```sh
#!/bin/sh

if [ -n "${xml_catalog_files_libxml2:-}" ]; then
    export XML_CATALOG_FILES="${xml_catalog_files_libxml2}"
else
    unset XML_CATALOG_FILES
fi
unset xml_catalog_files_libxml2
```
