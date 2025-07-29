#!/bin/bash


while getopts "d:t:r" opt
do 
    case $opt in 
        d)
            fastq_dir=$OPTARG
            ;;
        t)
            threshold=$OPTARG
            ;;
        r)
            remove=$OPTARG
            ;;
        \?)
            echo "Invalid option: $OPTARG"
            ;;
    esac
done


for i in "${fastq_dir}"/*fastq*
do 

    name=$(basename $i)
    size=$(du -sh $i | awk '{print $1}')
 
    if grep -Eq "[MGT]"  <<< $size
    then
        
        t_size=$( echo $size | awk '{gsub(/[^0-9]/, ""); print}')
        t_size=$(echo $size | awk '
            /G/ {printf "%.0f", $1 * 1024}
            /M/ {printf "%.0f", $1}
            /K/ {printf "%.0f", $1 / 1024}
            !/K|G|M/ {print $1}
        ')
       
       

        if [ $t_size -lt $threshold ]
        then

          
            if [ ! -d "${fastq_dir}"/small_files ]
            then
                mkdir "${fastq_dir}"/small_files
            fi
            mv $i "${fastq_dir}"/small_files
       
        fi
    
    elif grep -q "K" <<< $size
    then
        
        if [ ! -d "${fastq_dir}"/small_files ]
        then
            mkdir "${fastq_dir}"/small_files
        fi
        mv $i "${fastq_dir}"/small_files
    fi
done