
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
   git clone https://github.com/yourusername/bactflow.git

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

  





# Usage
## Run the Nextflow workflow directly from GitHub with the following command:

```sh
nextflow run bactflow --help -r main

Usage: nextflow run bactflow [options]

Options:
   
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
    --checkm_lineag_check   #If true, the genomes will be checked for their lineage completeness in one bin (default false).
    --genome_extension      #Required if '--checkm_lineag_check true'; default fasta.
    --run_flye              #If true, it runs Flye assembler; default true.
    --run_unicycler         #If true, it runs Unicycler hybrid assemlber, default false.
    --run_megahit           #If true, it runs Megahit assembler, default false.
    --run_spades            #If true, it runs Spades assembler, default false.
    --tax_class             #If true, it runs GTBtk taxonomic classification, default true.
    --prok_annot            #If true, it runs gene annotaiton by Prokka, default false. 
    --run_checkm            #If true, it runs checmk lineage and phylogenetic tree workflow.
    --checkm_db             #An absolute path to the Checkm database.  
    --gtdbtk_data_path      #Absolute path to the GDBtk database. 
    --run_quast             #Post-assembly stats by Quaset, default true.
    --genome_dir            #Path to already assembled genomes, only to run post-assembly tasks, e.g. taxonomy classification, gene annotations and quast or checkm

```

## Dependencies

<h3> <a href="https://github.com/Ecogenomics/GTDBTk/tree/master" target="_blank">GTDBtk</a></h3> 
You can run taxonomic classification workflow by setting `--tax_class true`. Alternatively, you can perform taxonomy classfiication on already-assembled genomes by adding `--run_flye false`. In this case you must provide the workflow with `--genome_dir`; a directory to the assembled genome and a genome `--extension fasta` or `fa`.  

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





