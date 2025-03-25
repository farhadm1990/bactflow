#!/bin/bash -ue
source $(conda info --base)/etc/profile.d/conda.sh




    if command -v mamba
    then
        conda_con=$(mamba env list | grep "bactflow" | awk '{print $1}' | grep -o '[a-zA-Z]' | wc -l)
        if [ $conda_con -eq 0 ]
        then
            mamba update --all -y
            mamba env create -f /home/hackion/Dropbox/Old-2gb/postdoc/ONT_helper/scripts/UI_bactflow/app/assem_app/config.yml
            bash /home/hackion/Dropbox/Old-2gb/postdoc/ONT_helper/scripts/UI_bactflow/app/assem_app/r_installer_pkg.sh
            echo "bactflow was successfully installed!"
            conda activate bactflow
        else 
            echo "Bactflow is already installed"
            conda activate bactflow

        fi
    else
        conda_con=$(conda env list | grep "bactflow" | awk '{print $1}' | grep -o '[a-zA-Z]' | wc -l)
        if [ $conda_con -eq 0 ]
        then
            conda update --all -y
            conda env create -f /home/hackion/Dropbox/Old-2gb/postdoc/ONT_helper/scripts/UI_bactflow/app/assem_app/config.yml
            bash /home/hackion/Dropbox/Old-2gb/postdoc/ONT_helper/scripts/UI_bactflow/app/assem_app/r_installer_pkg.sh
            echo "bactflow was successfully installed!"
            conda activate bactflow
        else 
            echo "Bactflow is already installed"
            conda activate bactflow
        fi
    fi



touch environment_created
