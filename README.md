<div style="text-align: center; margin-top: 0;">
    <img src="https://github.com/farhadm1990/bactflow/blob/main/logo/logo.png" alt="BactFlow Logo" width="400" height="400"/>
    <p><strong>BactFlow</strong> logo was designed by DALL·E :)</p>
</div>



# BactFlow

## Introduction

BactFlow is a workflow with a user-friendly UI for bacterial genome assembly of single-isolate sequencing reads from **Oxford Nanopore (ONT)**, **Illumina**, **PacBio**, and **ONT+Illumina hybrid** data. It is built with Nextflow DSL 2 and reads the usual outputs of Guppy and Dorado basecallers.

There are three modules:

| Module | What it does | UI port |
| --- | --- | --- |
| **Pre-assembly** | Concatenate barcode folders, drop tiny files, seqkit stats, read plots | 5000 |
| **Assembly** | Flye (ONT / PacBio), SPAdes (Illumina), Unicycler (hybrid), optional Medaka / circulate / QUAST / Bakta / GTDB-Tk / CheckM | 5002 |
| **Post-assembly** | QUAST, Bakta, GTDB-Tk, CheckM, SNP finder, Medaka variant calling, circular plots, strain finder | 5001 |

Only **one assembler** runs per assembly invocation. You can reuse the same `--out_dir` so Illumina, ONT, and PacBio FASTAs pile up in the same folder.

## Requirements

- [Docker](https://www.docker.com/) (**recommended** for the UI), **or**
- [Conda](https://docs.conda.io/en/latest/miniconda.html) / Mamba + [Nextflow](https://www.nextflow.io/docs/latest/index.html) for a local install

### Docker group (required once)

```sh
sudo groupadd docker
sudo usermod -aG docker $USER
newgrp docker
docker run hello-world
```

---

# Quick start (recommended)

This is the shortest path to open the UI successfully on a Linux desktop with Docker.

```sh
git clone https://github.com/farhadm1990/bactflow.git
cd bactflow

# Use an absolute work directory that can hold reads + outputs
WORK_DIR="$(pwd)"
mkdir -p "${WORK_DIR}/bactflow_out"

# Optional but useful for raw `docker run` browser pop-ups
./bactflow.sh --install-browser-hook

# Pre-assembly UI (port 5000) — browser opens when ready
./bactflow.sh preassem "${WORK_DIR}"

# In another terminal: assembly UI (port 5002)
./bactflow.sh assem "${WORK_DIR}" --cpus 10 --memory 16g

# Post-assembly UI (port 5001)
./bactflow.sh postassem "${WORK_DIR}" --cpus 10 --memory 16g
```

Point the UI forms at the bundled test reads under `${WORK_DIR}/test_reads/` (see below) and set **Output directory** to `${WORK_DIR}/bactflow_out` (or any absolute path **inside** `WORK_DIR`).

**Important:** Docker only sees paths under `WORK_DIR`. FASTQ folders, genome folders, databases (or symlinks to them), and `--out_dir` must all live under that absolute path.

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

| Module | Suggested test path (absolute) |
| --- | --- |
| Pre-assembly / Assembly ONT | `…/bactflow/test_reads/ont` |
| Assembly Illumina | `…/bactflow/test_reads/illumina` |
| Assembly PacBio | `…/bactflow/test_reads/pacbio` |
| Hybrid Unicycler | long: `…/test_reads/ont`, short: `…/test_reads/illumina` |

---

# Running the UI

## Docker: `bactflow.sh`

By default `bactflow.sh` **pulls the latest tag** from Docker Hub (`farhadm1990/bactflow_preassem`, `bactflow_assem`, `bactflow_postassem`). Use `--local` / `--rebuild` to build from this repo instead.

```sh
./bactflow.sh preassem /home/user/work_dir
./bactflow.sh assem /home/user/work_dir --cpus 16 --memory 32g
./bactflow.sh postassem /home/user/work_dir
```

Useful flags:

| Flag | Meaning |
| --- | --- |
| `--cpus N` / `--memory SIZE` | Resource limits (defaults: auto for preassem; 10 CPUs / 16g for assem & postassem) |
| `--port PORT` | Host port (defaults 5000 / 5002 / 5001) |
| `--no-browser` | Do not open a browser tab |
| `--install-browser-hook` | Install/start host helper so raw `docker run` can open the browser |
| `--local` | Use a local image build |
| `--rebuild` | Rebuild local image (implies `--local`) |
| `--tag TAG` | Pin a Hub tag (e.g. `v1.0`) |
| `--pull` | Force a Hub refresh |

```sh
./bactflow.sh assem /home/user/work_dir --cpus 16 --memory 32g --port 5002
./bactflow.sh preassem /home/user/work_dir --local --rebuild
```

The assembly image includes Flye, SPAdes, Unicycler, Circlator, and QUAST. Large taxonomy DBs and some polish/annotation tools are configured from the UI / post-assembly module.

### Raw `docker run` (optional)

Prefer `./bactflow.sh` when possible. If you run containers yourself:

1. Install the browser hook once: `./bactflow.sh --install-browser-hook`
2. Map a **fixed** host port (`-p 5000:5000`, not bare `-p 5000`)
3. Mount your home or work tree so the hook can see drop files if needed

```sh
docker run --rm --cpus=10 --memory=28g \
  --add-host=host.docker.internal:host-gateway \
  -p 5000:5000 \
  -v /home/user:/home \
  -v /home/user/work_dir:/home/user/work_dir \
  -e BACTFLOW_HOST_PORT=5000 \
  farhadm1990/bactflow_preassem:v1.0
```

<div style="text-align: center; margin-top: 10;">
    <img src="https://github.com/farhadm1990/bactflow/blob/main/pix/preassem.png" alt="preassem" style="max-width: 100%; height: auto;"/>
    <p><strong>Fig 1. </strong> Interface of the pre-assembly module. </p>
</div>
<br>

<div style="text-align: center; margin-top: 10;">
    <img src="https://github.com/farhadm1990/bactflow/blob/main/pix/assem.png" alt="Assem" style="max-width: 100%; height: auto;"/>
    <p><strong>Fig 2. </strong> Interface of the assembly module. </p>
</div>
<br>

<div style="text-align: center; margin-top: 10;">
    <img src="https://github.com/farhadm1990/bactflow/blob/main/pix/postassem.png" alt="Post assem" style="max-width: 100%; height: auto;"/>
    <p><strong>Fig 3. </strong> Interface of the post-assembly module. </p>
</div>

## Local UI (no Docker)

After the conda environment is installed (see [Installation](#installation-local-run-without-docker)):

```sh
conda activate bactflow
cd /path/to/bactflow

python3 UI/pre_assem_app/pre_assembly.py      # http://127.0.0.1:5000
python3 UI/assem_app/assembly.py              # http://127.0.0.1:5002
python3 UI/post_assem_app/post_assembly.py    # http://127.0.0.1:5001
```

Use absolute paths to `test_reads/…` and your output directory in the forms.

---

# Installation: local run without Docker

1. **Clone the repository and create the environment**

```sh
git clone https://github.com/farhadm1990/bactflow.git
cd bactflow
bash scripts/setup_bactflow_env.sh
```

After `conda activate bactflow`, Java and Nextflow are configured automatically (Java 17–24, Nextflow 24.10.2). If you have an older env, run:

```sh
mamba install -n bactflow openjdk=21 'nextflow=24.10.2'
bash scripts/install_bactflow_hooks.sh
conda deactivate && conda activate bactflow
nextflow -version
```

First-time package install only:

```sh
nextflow run UI/assem_app/main.nf --setup_only true --out_dir bactflow_out
```

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

All examples assume:

```sh
conda activate bactflow
cd /path/to/bactflow
```

Help for the current assembly workflow:

```sh
nextflow run UI/assem_app/main.nf --help
```

Add `-resume` to continue an interrupted run. Use `--cpus` for threads. Set `--out_dir` to an absolute path if you can.

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

From the post-assembly web UI (port **5001**) you can also:

- **SNP finder** — BWA + bcftools; interactive VCF tables per genome under `out_dir/snps/`
- **Variant calling (Medaka)** — tables from each genome’s **`medaka.vcf`** under `out_dir/vcs/`
- **Circular plot** — from Bakta GBK/GFF outputs
- **Strain finder** — abundance / prevalence tables and plots from enzyme + annotation inputs

---

# Output directory (`bactflow_out`)

Published results land here. Nextflow work files go to `bactflow_out/.nextflow-work` and are removed after a successful assembly run.

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

- Prefer `./bactflow.sh …` (opens the host browser after the UI is ready).
- For raw `docker run`, run `./bactflow.sh --install-browser-hook` once and use `-p HOST:CONTAINER` (e.g. `-p 5000:5000`).
- Rebuild local images after UI changes: `./bactflow.sh <module> <work_dir> --local --rebuild`.

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
