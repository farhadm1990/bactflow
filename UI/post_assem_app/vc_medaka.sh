#!/bin/bash



display_help(){
        echo "Usage: $0 -g <folder containing genomes> -r <reference genome> -o <output directory> -c <number of cpus> "
        echo "  -g  Genomes directory (folder containing the genomes to be compared to the reference.)"
        echo "  -r  Reference genome to for variants to be called against. .fasta or .fa"
        echo "  -o  Output directory."
        echo "  -c  Number of cpus." 

}

while getopts ":hg:r:o:c:" opt #hg, is not seperated to tell that h doesn't need an argument
do
    case $opt in 
        h)
            display_help
            exit 0
            ;;

        g)
            genomes="$OPTARG"
            ;;
        r)
            reference="$OPTARG"
            ;;
        o)
            output="$OPTARG"
            ;;
        c)
            ncpu="$OPTARG"
            ;;
        \?)
            echo "Inavlid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument! " >&2
            exit 1
            ;;
    esac
done

if [ "$#" -eq 0 ] || [ -z "$genomes" ] || [ -z "$reference" ] || [ -z "$output" ] || [ -z "$ncpu" ]
then 
    echo "Missing required argument(s)! " >&2
    display_help
    exit 1
fi

if [ ! -d "$output" ]
then 
    mkdir -p "$output"
fi 


# source activate medaka
for genome in "$genomes"/*.fasta
do 
    if [ ! -d "${output}" ]
    then
        mkdir -p "${output}"
    fi

    name=$(basename $genome | cut -f 1 -d'.' )
    ref_name="${reference}"
    refname=$(echo $ref_name | cut -f1 -d'.')
    refname="$(echo ${refname##*_})"
    out_name=$(echo $name"_vs_"$refname)
    medaka_variant -t $(($ncpu)) -o "${output}"/$out_name -i $genome -r $ref_name
    rm -rf *.mmi *.fai
done
