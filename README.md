# Genebrosh

## Introduction

Genebrosh is a workflow for bacterial genome assembly of single isolate and metagenomics sequencing reads extracted from Oxford Nanopore Technology (ONT) and Illumina platforms. It is designed using Nextflow DSL 2 technology and reads the generic outputs of Guppy and Dorado basecallers.

## Requirements

- [Nextflow](https://www.nextflow.io/docs/latest/index.html)
- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Docker](https://www.docker.com/)

## Installation

1. **Clone the Repository**
   ```sh
   git clone https://github.com/yourusername/genebrosh.git
   cd genebrosh
Set Up Conda Environment
sh
Copy code
conda create -n ont_helper -f environment.yml
conda activate ont_helper
````


# Usage
## Run the Nextflow workflow directly from GitHub with the following command:

```sh
nextflow run yourusername/genebrosh \
  --fastq_dir /path/to/fastq \
  --extension '.fastq.gz' \
  --cpus 30 \
  --tensor_batch 200 \
  --genome_size 6 \
  --coverage 80 \
  --tax_class true \
  --checkm_lineage_check true \
  --prok_annot true

  ```

