#!/bin/bash
# set -x 
set -e 

display_help(){
    echo "...::: Variant finder :::..."
    echo "Usage $0 -r <reference fasta file> -g <genome | query or a directory with genomes in fasta> -t <n cpus> -o <output dir>"
    echo "  -r  Reference genome to for variants to be called against. .fasta or .fa"
    echo "  -g  Genomes directory (folder containing the genomes to be compared to the reference.)"
    echo "  -o  Output directory."
    echo "  -t  Number of cpus." 
}


while getopts ":hg:r:g:o:t:" opt
do 
    case $opt in 
        h)
            display_help
            exit 0
            ;;
        r)
            reference="$OPTARG"
            ;;
        g)
            genomes="$OPTARG"
            ;;
        o)
            output_dir="$OPTARG"
            ;;
        t)
            cpus="$OPTARG"
            ;;
        \?)
            echo "Invalid option. Please try again"
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument" >&2
            exit 1
            ;;


    esac
done


if [ "$#" -eq 0 ] || [ -z "$genomes" ] || [ -z "$reference" ] || [ -z "$output_dir" ] || [ -z "$cpus" ]
then 
    echo "Missing required argument(s)! " >&2
    display_help
    exit 1
fi

genomes=$(realpath "$genomes")
output_dir=$(realpath "$output_dir")

#export "$genomes" 
#export "$output_dir" 
#export "$reference"
#export "$cpus"

if [ ! -d "$output_dir" ]
then 
    mkdir -p "$output_dir"
fi 

#index reference
echo "Indexing the reference $reference"

bwa index "${reference}"

for genome in "${genomes}"/*
do 
    name=$(basename "$genome" | cut -f 1 -d'.' )
    if [ ! -d "${output_dir}/${name}" ]
    then
        mkdir -p "${output_dir}/${name}"
    fi

    
    echo "Processing $genome ...."

    refname=$(echo "${reference}" | cut -f1 -d'.')
    refname="$(echo ${refname##*_})"
    out_name=$(echo $name"_vs_"$refname)
    bwa mem "${reference}" $genome -t $(($cpus)) -o "${output_dir}/${name}/${out_name}.sam"
    
    # sam to bam

    samtools view -Sb "${output_dir}/${name}/${out_name}.sam" > "${output_dir}/${name}/${out_name}.bam"
    rm -f "${output_dir}/${name}/${out_name}.sam"

    samtools sort "${output_dir}/${name}/${out_name}.bam" > "${output_dir}/${name}/${out_name}"_al.bam
    samtools index "${output_dir}/${name}/${out_name}"_al.bam

    # calling variant

    bcftools mpileup -Ou -f "${reference}" "${output_dir}/${name}/${out_name}"_al.bam --threads $(($cpus)) | bcftools call -mv -Ov -o "${output_dir}/${name}/${out_name}".vcf

    echo "Variant calls for $genome is ready at '${output_dir}/${name}'"
done


rm -rf *.mmi *.fai *.bai *.amb *.pac *.sa



#aling genome to the ref






