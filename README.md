
<div style="text-align: center; margin-top: 0;">
    <img src="https://github.com/farhadm1990/bactflow/blob/main/logo/logo.png" alt="BactFlow Logo" width="400" height="400"/>
    <p><strong>BactFlow</strong> logo was desinged by DALLE :) </p>
</div>



# BactFlow

## Introduction

BactFlow is a workflow for bacterial genome assembly of single isolate and metagenomics sequencing reads extracted from Oxford Nanopore Technology (ONT) and Illumina platforms. It is designed using Nextflow DSL 2 technology and reads the generic outputs of Guppy and Dorado basecallers.

## Requirements

- [Nextflow](https://www.nextflow.io/docs/latest/index.html)
- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Docker](https://www.docker.com/)

## Installation

1. **Clone the Repository**
```sh
   git clone https://github.com/farhadm1990/bactflow.git

   cd bactflow
   mamba env create -f config.yml

   cd ..
   cp bactflow ~/.nextflow/assets/farhadm1990/

   nextflow run bactflow --help
```

2. **Clone the Repository via `nextflow pull` command**
```sh

nextflow pull farhadm1990/bactflow

nextflow  run bactflow
```

3. **Setting up the conda environment by `bactflow` internal funciton `envSetUp()`:**
On the first lunch you can run the following code to create conda environment called bactflow and install the required packages within it.
```sh
nextflow run bactflow -r main --setup_only true
```
NOTE: the option `--setup_only` is by default `false` which means, the `envSetUp()` function will be invoced automatically before running the downstream processes. 
  


# Usage
## Run the Nextflow workflow directly from GitHub with the following command:

```sh
nextflow run bactflow --help -r main

Usage: nextflow run bactflow [options]

Options:
    --setup_only            #Only runs envSetUp(), default false.
    --fastq_dir             #Absolute path to the fastq_pass directory (required). 
    --concat_reads          #Default true, it concatenates all your ONT basecaller 4000-chunk reads into one fastq file. Set it to false if it is already concatenated.
    --extension             #String; extention of basecalled fastq files; default '.fastq.gz'
    --cpus                  #Number of available cpus; default 1.
    --coverage_filter       #If you want to normalize all your genomes to a certain coverage (default false).
    --coverage              #Only if '--coverage_filter true'; default is 50.
    --genome_size           #Genome size for coverage normalizaiton. Only if '--coverage_filter true'; default is 6.
    --out_dir               #Output directory of your final results. Default "genebrosh_output"
    --tensor_batch          #Medaka tensorflow batch size. Lower it in low coverage genomes. Default 200.
    --nanofilter            #Filtering reads for length and quality; default true.
    --min_length            #If '--nanofilter' true, filter reads below a certain read length (default 1000). 
    --min_quality           #If '--nanofilter' true, filter reads below a certain read quality (default 16 for R10.4.1 flowcells). 
    --medaka_polish         #If true, it will polish assembled genomes by medaka (dfault false).
    --basecaller_model      #Basecaller model for medaka polishing step. 'r1041_e82_400bps_hac_v4.2.0'
    --genome_extension      #Required if '--checkm_lineag_check true'; default fasta.
    --run_flye              #If true, it runs Flye assembler; default true.
    --circle_genome         #If true, it will fix the start of the genome to an arbitrary gene, e.g. dnaA. Default false.
    --run_unicycler         #If true, it runs Unicycler hybrid assemlber, default false.
    --run_megahit           #If true, it runs Megahit assembler, default false.
    --run_spades            #If true, it runs Spades assembler, default false.
    --tax_class             #If true, it runs GTBtk taxonomic classification, default false.
    --prok_annot            #If true, it runs gene annotaiton by Prokka, default false. 
    --run_checkm            #If true, it runs checmk lineage and phylogenetic tree workflow.
    --checkm_db             #An absolute path to the Checkm database.  
    --gtdbtk_data_path      #Absolute path to the GDBtk database. 
    --run_quast             #Post-assembly stats by Quaset, default true.
    --genome_dir            #Path to already assembled genomes, only to run post-assembly tasks, e.g. taxonomy classification, gene annotations and quast or checkm

```

## Dependencies

<h3> <a href="https://github.com/Ecogenomics/GTDBTk/tree/master" target="_blank">GTDBtk</a></h3> 
<p>You can run taxonomic classification workflow by setting <code>--tax_class true</code>. Alternatively, you can perform taxonomy classfiication on already-assembled genomes by adding <cod>--run_flye false</code>. In this case you must provide the workflow with <code>--genome_dir</code>; a directory to the assembled genome and a genome <code>--extension fasta</code> or <code>fa</code>. GTDBtk depends on an external <a href="https://ecogenomics.github.io/GTDBTk/installing/index.html">datasbase</a> which can be downloaded and extracted as follows:  </p>

```sh
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz

tar xvzf gtdbtk_data.tar.gz
```

<p>The full path of the extracted must be proveded to the pipline <code>nextflow.config</code> file or when running the pipeline via <code>--gtdbtk_data_path /home/databases/gtdbtk_db/release220</code> for example. 
</p>

# Generic ouputs 

## Report maker

<div style="text-align: center; margin-top: 10;">
    <img src="https://github.com/farhadm1990/bactflow/blob/main/pix/pan_genome.png" alt="Pangenome" style="max-width: 100%; height: auto;"/>
    <p><strong>Fig 1. </strong> Pangenome of assembled genomes extracted from bactflow. </p>
</div><br>

<div style="text-align: center; margin-top: 10;">
    <img src="https://github.com/farhadm1990/bactflow/blob/main/pix/tree_genome.png" alt="Phylogenetic tree" style="max-width: 100%; height: auto;"/>
    <p><strong>Fig 2. </strong> Phylogenetic tree of assembled genomes based on their ANI. </p>
</div>


## Potential issues

In the first lunch of the program, process `envSetUP` will be invoked and it creates a conda environment called `bactflow`. In some conda environment, the following error could be raised:
```sh
miniconda3/envs/bactflow/etc/conda/deactivate.d/libxml2_deactivate.sh: line 3: xml_catalog_files_libxml2: unbound variable
```
### Solution 
You can open the file and add edit it as follows and rerun `bactflow`:
```sh
#!/bin/sh

if [ -n "${xml_catalog_files_libxml2:-}" ]; then
    export XML_CATALOG_FILES="${xml_catalog_files_libxml2}"
else
    unset XML_CATALOG_FILES
fi
unset xml_catalog_files_libxml2


```

