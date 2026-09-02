
<div style="text-align: center; margin-top: 0;">
    <img src="https://github.com/farhadm1990/bactflow/blob/main/logo/logo.png" alt="BactFlow Logo" width="400" height="400"/>
    <p><strong>BactFlow</strong> logo was designed by DALL·E :)</p>
</div>



# BactFlow

## Introduction

BactFlow is a workflow with a user-friendly UI for bacterial genome assembly of single-isolate sequencing reads from Oxford Nanopore (ONT), Illumina, PacBio, and ONT+Illumina hybrid data. It is built with Nextflow DSL 2 and reads the usual outputs of Guppy and Dorado basecallers.

There are three modules:

| Module | What it does | UI port |
| --- | --- | --- |
| **Pre-assembly** | Concatenate barcode folders, drop tiny files, seqkit stats, read plots | 5000 |
| **Assembly** | Flye (ONT), SPAdes (Illumina), Unicycler (hybrid), Flye (PacBio), optional circulate / QUAST / Bakta / GTDB-Tk / CheckM | 5002 |
| **Post-assembly** | QUAST, Bakta, GTDB-Tk, CheckM on already-assembled FASTA files | 5001 |

Only **one assembler** runs per assembly invocation. You can reuse the same `--out_dir` so Illumina, ONT, and PacBio FASTAs pile up in the same folder.

## Requirements

- [Nextflow](https://www.nextflow.io/docs/latest/index.html)
- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Docker](https://www.docker.com/)

### NOTE: if you want to run the UI through Docker (recommended), remember to have Docker installed and your user in the `docker` group

```sh
sudo groupadd docker
sudo usermod -aG docker $USER
newgrp docker
# test
docker run hello-world
```

---

# Running the UI

## Docker (recommended): `bactflow.sh`

Provide a module name (`preassem`, `assem`, or `postassem`) and an **absolute** path to your working directory. Docker can only read and write inside that path, so your FASTQ folders, genome folders, databases (or symlinks to them), and `--out_dir` must all live under it.

```sh
./bactflow.sh preassem /home/user/work_dir
./bactflow.sh assem /home/user/work_dir
./bactflow.sh postassem /home/user/work_dir
```

Optional flags: `--cpus`, `--memory`, `--port`, `--no-browser`.

```sh
./bactflow.sh assem /home/user/work_dir --cpus 16 --memory 32g --port 5002
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

python3 UI/pre_assem_app/pre_assembly.py      # http://127.0.0.1:5000
python3 UI/assem_app/assembly.py              # http://127.0.0.1:5002
python3 UI/post_assem_app/post_assembly.py    # http://127.0.0.1:5001
```

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

Already one FASTQ per sample (`--concat_reads false`):

```text
ont_reads/
├── TL110.fastq.gz
├── TL19.fastq.gz
└── TL29.fastq.gz
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
├── TL110_R1.fastq.gz
├── TL110_R2.fastq.gz
├── TL19_Illumina_R1.fastq.gz
└── TL19_Illumina_R2.fastq.gz
```

For hybrid Unicycler, the ONT filename is peeled down to a sample id (`TL110`) and matched to Illumina files such as `TL110_*_R1.fastq.gz`.

### PacBio (Flye)

```text
pacbio/
├── TL110.fastq.gz
└── TL19.fastq.gz
```

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

### ONT — Flye

Raw ONT (default `--ont_read_type nano-raw`). Use `nano-hq` for Q20+/Dorado HAC/SUP, or `nano-corr` for already-corrected reads.

```sh
nextflow run UI/assem_app/main.nf \
  --run_flye true \
  --run_spades false \
  --run_unicycler false \
  --run_pacbio false \
  --fastq_dir /home/user/work_dir/ont_reads \
  --concat_reads false \
  --ont_read_type nano-raw \
  --nanofilter true \
  --min_length 1000 \
  --min_quality 16 \
  --circle_genome true \
  --run_quast true \
  --cpus 10 \
  --out_dir bactflow_out \
  -resume
```

Optional Medaka polish (needs the basecaller model that matches the reads):

```sh
nextflow run UI/assem_app/main.nf \
  --run_flye true \
  --fastq_dir /home/user/work_dir/ont_reads \
  --concat_reads false \
  --ont_read_type nano-hq \
  --medaka_polish true \
  --basecaller_model r1041_e82_400bps_hac_v4.2.0 \
  --circle_genome true \
  --run_quast true \
  --cpus 10 \
  --out_dir bactflow_out
```

Coverage downsample before Flye:

```sh
nextflow run UI/assem_app/main.nf \
  --run_flye true \
  --fastq_dir /home/user/work_dir/ont_reads \
  --concat_reads false \
  --coverage_filter true \
  --coverage 50 \
  --genome_size 6 \
  --out_dir bactflow_out \
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
  --fastq_dir /home/user/work_dir/illumina \
  --circle_genome true \
  --run_quast true \
  --cpus 10 \
  --out_dir bactflow_out \
  -resume
```

### Hybrid ONT + Illumina — Unicycler

Long reads in `--fastq_dir`, short-read pairs in `--short_read_dir`. Sample prefixes must match (for example ONT `TL110_Native_R1041.fastq` with Illumina `TL110_Illumina_R1.fastq.gz`).

```sh
nextflow run UI/assem_app/main.nf \
  --run_flye false \
  --run_spades false \
  --run_unicycler true \
  --run_pacbio false \
  --fastq_dir /home/user/work_dir/ont_reads \
  --short_read_dir /home/user/work_dir/illumina \
  --concat_reads false \
  --circle_genome true \
  --run_quast true \
  --cpus 10 \
  --out_dir bactflow_out \
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
  --fastq_dir /home/user/work_dir/pacbio \
  --concat_reads false \
  --pacbio_read_type pacbio-hifi \
  --circle_genome true \
  --run_quast true \
  --cpus 10 \
  --out_dir bactflow_out \
  -resume
```

### Assembly plus taxonomy / annotation / CheckM

Point databases at absolute paths (or at symlinks under the Docker work directory).

```sh
nextflow run UI/assem_app/main.nf \
  --run_flye true \
  --fastq_dir /home/user/work_dir/ont_reads \
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
  --out_dir bactflow_out
```

## 3. Post-assembly (already have FASTA files)

Skip assemblers and score genomes in `--genome_dir` (usually `bactflow_out/asm_out_dir/fastas` or `bactflow_out/circulated_fasta`).

```sh
nextflow run UI/post_assem_app/main.nf \
  --run_flye false \
  --genome_dir /home/user/work_dir/bactflow_out/asm_out_dir/fastas \
  --genome_extension fasta \
  --out_dir bactflow_out \
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
└── checkm_out/                               # --run_checkm true (optional)
    ├── checkm_lineage.txt
    ├── taxon_tree.newick
    ├── genome_tree.newick
    └── genome_tree.tree
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

On the first launch, process `envSetUP` creates a conda environment called `bactflow`. In some conda setups this error can appear:

```sh
miniconda3/envs/bactflow/etc/conda/deactivate.d/libxml2_deactivate.sh: line 3: xml_catalog_files_libxml2: unbound variable
```

### Solution

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
