#!/bin/bash

source $(conda info --base)/etc/profile.d/conda.sh

conda activate ont_helper


R -e "BiocManager::install("phyloseq")"
R -e "BiocManager::install("ggtree")"
R -e "install.packages("patchwork")"
R -e "install.packages("readr")"
R -e "install.packages("ggpubr")"