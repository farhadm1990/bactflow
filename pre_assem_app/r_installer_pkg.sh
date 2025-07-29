#!/bin/bash

source $(conda info --base)/etc/profile.d/conda.sh

conda activate bactflow
# sudo apt update && sudo apt install -y libcurl4-openssl-dev  pkg-config

R -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); \
if(!require('phyloseq')) {BiocManager::install('phyloseq',  ask=FALSE)}; \
if(!require('Biostrings')) {BiocManager::install('Biostrings',  ask=FALSE)}; \
if(!require('seqinr')) {install.packages('seqinr')}; \
if(!require('msa')) {BiocManager::install('msa', force = TRUE,  ask=FALSE)}; \
if(!require('ape')) {install.packages('ape')}; \
if(!require('tidyverse')) {install.packages('tidyverse')}; \
if(!require('patchwork')) {install.packages('patchwork')}; \
if(!require('readr')) {install.packages('readr')}"



if [ $? -eq 0 ]
then 
    echo "All R packages have been installed!"

fi