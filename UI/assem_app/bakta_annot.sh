#!/bin/bash

genomes="./"
cpus=35
out_dir="bakta_out"
db_dir=""

display_help(){
    echo "Usage: $0 -g <genomes_directory> -c <cpus> -o <output_directory>" #$0 is the name of the script
    echo "Options:"
    echo "  -g   Genomes directory to fasta files (default: ./)"
    echo "  -c   Number of CPUs (default: 35)"
    echo "  -d   Directory to bakta database"
    echo "  -o   Output directory (default: out_dir)"
    exit 0
}

while getopts ":g:c:d:o:" opt
do 
    case $opt in
        g) 
            genomes="$OPTARG"
            ;;
        c)
            cpus="$OPTARG"
            ;;
        d)
            db_dir="$OPTARG"
            ;;
    
        o) 
            out_dir="$OPTARG"
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

if [ "$#" -eq 0 ]
then 
    display_help
fi

if [ ! -d "$out_dir" ]
then 
    mkdir -p "$out_dir"
fi 




cpus=$(($cpus))
for file in $(ls "$genomes"/*.fasta)
do 
    striped="$(basename $file | cut -d'.' -f1)"
    

    dest="${out_dir}"/"$striped"_bakta
  
 
    bakta --db "${db_dir}" "${file}" --prefix "$striped" --output "$dest" --force --threads "${cpus}"
done

