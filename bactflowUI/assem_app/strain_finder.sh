#!/bin/bash



set -e 


display_help(){
    echo "Usage: $0 quantifies annotated genes in different genomes"
    echo ""
    echo "-d <directory_gene_files> -g <gene_type> -f <out_file_name> -p <prevalance: default true> -o <output_directory> -n <name_of_organism> -e <enzyme_table>"
    echo
    echo "Options:"
    echo "  -d   Directory with gene annotation tab-delimited files"
    echo "  -g   Type of genes: default cds"
    echo "  -f   Output name for gene count table: default gene_count_table"
    echo "  -c   Only create gene count table: default true"
    echo "  -p   Prevalence or abundance of genes? (default: true i.e. prevalence)"
    echo "  -o   Output directory (default: gene_finder_out)"
    echo "  -n   Name of the organism, e.g. L. plantarum (default: Bacteria)"
    echo "  -e   Enzyme file: a tab-delimited file with enzymes of your interest"
    echo "  -w   Plot width (default: 20)"
    echo "  -l   Plot length (default: 20)"
    echo "  -h   Print this help message and exit"
    exit 0
}

prevalance="true"
gene_type="cds"
out_name="gene_count_table"
out_dir="gene_finder_out"
organism="Bacteria"
width=$((20))
height=$((20))
count_only='true'

while getopts ":d:g:f:c:p:o:n:e:w:l:h" args
do  
    case $args in
        d)
            gene_dir="$OPTARG"
            ;;

        g) 
            gene_type="$OPTARG"
            ;;
        f)
            out_name="$OPTARG"
            ;;
        c)
            count_only="$OPTARG"
            ;;
        p)
            prevalance="$OPTARG"
            ;;
        
        o)
            out_dir="$OPTARG"
            ;;
        
        n)
            organism="$OPTARG"
            ;;

        e)
            enzyme_file="$OPTARG"
            ;;
        w)
            width="$OPTARG"
            ;;
        l)
            height="$OPTARG"
            ;;

        h)

            display_help
            exit 0
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



chmod +x ./gene_counter_bakta.py
chmod +x ./gene_plotter.r

if [ $prevalance = "true" ]
then 
    prevalance="TRUE"
elif [ $prevalance = "false" ]
then
    prevalance="FALSE"
fi 

script_dir="$(dirname "$(realpath "$0")")"

./gene_counter_bakta.py -d "${gene_dir}" -t "${gene_type}" -o "${out_dir}/${out_name}"

if [ $count_only = 'false' ]
then 

    in_name=$(echo "${out_dir}/${out_name}.tsv")

   ./gene_plotter.r -c "${in_name}" -p "${prevalance}" -o "${out_dir}" -n "${organism}" -e "${enzyme_file}" -w "${width}" -l "${height}"

fi