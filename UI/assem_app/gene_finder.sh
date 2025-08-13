#!/bin/bash

help(){
    echo "Usage $0: finds the position of genes on chromosomes::.... "
    echo 
    echo "-d    Directory of tab-seperated gene files, tsv"
    echo "-g    List of comma-seperated (without space) genes to be located"
    echo "-o    Output directory; default curent dir"
}

genes=""

while getopts ":d:g:o:h" opt
do  
    case $opt in
    d)
        gene_dir="$OPTARG"
        ;;
    g)
        genes="$OPTARG"
        ;;
    o)
        out_dir="$OPTARG"
        ;;
    h)
        help
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


colname="Gene\tstart_position\tend_position\tStrand"

if [ "$#" -eq 0 ]
then    
    help
fi


if [ ! -d "$out_dir" ]
then 
    mkdir -p "$out_dir"
    echo "$out_dir has been created!"
else
    echo "$out_dir already exists and I proceed to the next step!"
fi 

files=$(find $gene_dir -type f -name "*.tsv")

for fold in $files
do
    if [ -f $fold ]
    then 
        name=$(basename $fold '.tsv')
        echo -e "$colname" > "${out_dir}"/$name.tsv

        grep -E "$(echo $genes | sed 's/,/|/g')"  $fold | awk '{print $7 "\t" $3 "\t" $4 "\t" $5}' >> "${out_dir}"/$name.tsv
    else
        echo "Error: there is not gene TSV file in this directory. Please try again!"
        exit 1
    fi
done