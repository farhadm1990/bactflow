#!/bin/bash

source $(conda info --base)/etc/profile.d/conda.sh

conda activate bactflow


R -e "BiocManager::install('phyloseq')"
R -e "BiocManager::install('Biostrings')"
R -e "BiocManager::install('ggtree')"
R -e "BiocManager::install('msa')"
R -e "install.packages('ape')"
R -e "BiocManager::install('tidyverse')"
R -e "install.packages('patchwork')"
R -e "install.packages('readr')"
R -e "install.packages('ggpubr')"