#!/bin/bash

source $(conda info --base)/etc/profile.d/conda.sh

conda activate bactflow
# sudo apt update && sudo apt install -y libcurl4-openssl-dev  pkg-config

R -e "options(repos = c(CRAN = 'https://cloud.r-project.org'))"
R -e "if(!require('phyloseq')) {BiocManager::install('phyloseq')}"
R -e "if(!require('Biostrings')) {BiocManager::install('Biostrings')}"
#R -e "if(!require('ggtree')) {BiocManager::install('ggtree')}"
R -e "if(!require('seqinr')) {install.packages('seqinr')}"
R -e "if(!require('msa')) {BiocManager::install('msa', force = TRUE)}"
R -e "if(!require('ape')) {install.packages('ape')}"
R -e "if(!require('tidyverse')) {install.packages('tidyverse')}"
R -e "if(!require('patchwork')) {install.packages('patchwork')}"
R -e "if(!require('readr')) {install.packages('readr')}"