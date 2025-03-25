#!/bin/bash -ue
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bactflow
conda env list | grep bactflow > checked.txt
