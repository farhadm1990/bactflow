
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
   cp bactflow ~/.nextflow

   nextflow run bactflow
```
  





# Usage
## Run the Nextflow workflow directly from GitHub with the following command:

   ```sh
   nextflow run farhadm1990/bactflow \ 
   --fastq_dir /path/to/fastq \         # Full path to the fastq_pass dir
   --extension '.fastq.gz' \            # extension of fastq reads 
   --cpus 30 \                          # threads
   --tensor_batch 200 \                 # number of tensorflow batch for medaka polishing. Consider reducing it in low-coverage samples.
   --genome_size 6 \                    
   --coverage 80 \
   --tax_class true \
   --checkm_lineage_check true \
   --prok_annot true \
   --run_quast true \
   --medaka_polish true



