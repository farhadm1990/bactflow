#!/bin/bash

extension=".fastq.gz"
cpus=1
output_dir="bactflow_out"

while getopts :g:c:e:o: opt
do  
    case $opt in
        g) 
            fastq_dir="$OPTARG" 
            ;;
        c)
            cpus="$OPTARG"
            ;;
        e)
            extension="$OPTARG"
            ;;
        o)
            output_dir="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done


source $(conda info --base)/etc/profile.d/conda.sh
conda activate bactflow

if [ "$#" -eq 0 ]
then 
    echo "Usage: $0 -g <fastq_directory> -c <cpus> -e <extension>"
    echo "Options:"
    echo "  -g   Fastq directory (default: ./)"
    echo "  -c   Number of CPUs (default: 1)"
    echo "  -e   Extension of the fastq files (default: .fastq.gz)"
    echo "  -o   Output directory (default: ./)"
    exit 0
fi

if [ -z "${fastq_dir}" ]
then
    echo "no directory is specified to $fastq_dir"
    eixt 1
fi
   fastq_dir=$(realpath "${fastq_dir}") 

    num_file=$(ls "${fastq_dir}"| wc -l)

    if [ -d "${fastq_dir}"/pooled ] 
    then 
        num_pooled=$(ls "${fastq_dir}"/pooled | wc -l)
        num_dif=$(("${num_file}" - "${num_pooled}"))
        if [ $num_dif -eq 1 ]
        then
           touch concatenated_fq_are_ready
        else 
            if [ "${extension}" == ".fastq.gz" ]
            then
                
                parallel --will-cite -j   "${cpus}" '
                    if [ -d {} ]  && [ "$(basename {})" != "pooled" ]
                    then
                        fastq_dir=$(dirname {})
                        name=$(basename {})
                       
                        zcat {}/*.fastq.gz >> "${fastq_dir}"/"${name}_pooled.fastq"
                        mv "${fastq_dir}"/"${name}"_pooled.fastq "${fastq_dir}/pooled"
                    fi
                '  ::: "${fastq_dir}"/*

                touch concatenated_fq_are_ready
                exit 0

            elif [ "${extension}" == ".fastq" ] 
            then	
              
                parallel --will-cite -j   "${cpus}" '
                    if [ -d {} ]  && [ "$(basename {})" != "pooled" ]
                    then
                        fastq_dir=$(dirname {})
                        name=$(basename {})
                        cat {}/*.fastq >> "${fastq_dir}"/"${name}_pooled.fastq"
                        mv "${fastq_dir}"/"${name}"_pooled.fastq "${fastq_dir}/pooled"
                    fi
                ' ::: "${fastq_dir}"/*
                
                touch concatenated_fq_are_ready
                exit 0
            else
                echo "Your extention is not recognized!"
                exit 1

            fi
        fi



    else 
        
        mkdir -p "${fastq_dir}"/pooled

        if [ "${extension}" == ".fastq.gz" ]
        then

          
            parallel --will-cite -j   "${cpus}" '
                if [ -d {} ]  && [ "$(basename {})" != "pooled" ]
                then
                    fastq_dir=$(dirname {})
                    name=$(basename {})
                    zcat {}/*.fastq.gz >> "${fastq_dir}"/"${name}_pooled.fastq"
                    mv "${fastq_dir}"/"${name}"_pooled.fastq "${fastq_dir}/pooled"
                fi
            ' ::: "${fastq_dir}"/*

            touch concatenated_fq_are_ready

        elif [ "${extension}" == ".fastq" ] 
        then

            export "${fastq_dir}"
            parallel --will-cite -j   "${cpus}" '
                if [ -d {} ]  && [ "$(basename {})" != "pooled" ]
                then
                    fastq_dir=$(dirname {})
                    name=$(basename {})
                    cat {}/*.fastq >> "${fastq_dir}"/"${name}_pooled.fastq"
                    mv "${fastq_dir}"/"${name}"_pooled.fastq "${fastq_dir}/pooled"
                fi
            ' ::: "${fastq_dir}"/*	
            
            touch concatenated_fq_are_ready
        else
            echo "Your extention is not recognized!"
            exit 1

        fi

    fi