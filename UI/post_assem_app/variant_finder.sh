#!/bin/bash
set -euo pipefail

display_help(){
    echo "...::: Variant finder :::..."
    echo "Usage $0 -r <reference fasta file> -g <directory with genomes in fasta> -t <n cpus> -o <output dir>"
    echo "  -r  Reference genome for variants to be called against (.fasta / .fa / .fna)"
    echo "  -g  Genomes directory (folder containing the genomes to compare to the reference)"
    echo "  -o  Output directory"
    echo "  -t  Number of cpus"
}


while getopts ":hg:r:o:t:" opt
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


if [ "${#}" -eq 0 ] || [ -z "${genomes:-}" ] || [ -z "${reference:-}" ] || [ -z "${output_dir:-}" ] || [ -z "${cpus:-}" ]
then
    echo "Missing required argument(s)!" >&2
    display_help
    exit 1
fi

if [ ! -f "$reference" ]
then
    echo "Reference file not found: $reference" >&2
    exit 1
fi

if [ ! -d "$genomes" ]
then
    echo "Genomes directory not found: $genomes" >&2
    exit 1
fi

reference=$(realpath "$reference")
genomes=$(realpath "$genomes")
output_dir=$(realpath "$output_dir")
mkdir -p "$output_dir"

is_fasta() {
    local path="$1"
    local base low
    base=$(basename "$path")
    low=$(printf '%s' "$base" | tr '[:upper:]' '[:lower:]')
    case "$low" in
        *.fasta|*.fa|*.fna) return 0 ;;
        *) return 1 ;;
    esac
}

bwa_index_complete() {
    local ref="$1"
    local ext
    for ext in amb ann bwt pac sa; do
        if [ ! -s "${ref}.${ext}" ]; then
            return 1
        fi
    done
    return 0
}

ensure_reference_index() {
    local ref="$1"
    # Keep BWA/samtools indexes under the SNP out dir — never beside the user's
    # original FASTA (avoids races and "fail to locate the index files").
    if bwa_index_complete "$ref" && [ -s "${ref}.fai" ]; then
        return 0
    fi
    echo "Indexing reference for BWA/bcftools: $ref"
    rm -f "${ref}.amb" "${ref}.ann" "${ref}.bwt" "${ref}.pac" "${ref}.sa" "${ref}.fai"
    bwa index "$ref"
    samtools faidx "$ref"
    if ! bwa_index_complete "$ref"; then
        echo "Failed to build BWA index for $ref" >&2
        exit 1
    fi
}

# Collect only real FASTA genomes; skip BWA/samtools index sidecars and the reference itself.
mapfile -t genome_files < <(
    find "$genomes" -maxdepth 1 -type f \( -iname '*.fasta' -o -iname '*.fa' -o -iname '*.fna' \) | sort
)

if [ "${#genome_files[@]}" -eq 0 ]
then
    echo "No .fasta/.fa/.fna files found in ${genomes}" >&2
    exit 1
fi

# Local working copy/symlink of the reference inside the output tree.
ref_work="${output_dir}/.ref_index"
mkdir -p "$ref_work"
ref_local="${ref_work}/reference.fasta"
ln -sfn "$reference" "$ref_local"
ensure_reference_index "$ref_local"

processed=0
for genome in "${genome_files[@]}"
do
    genome=$(realpath "$genome")
    if [ "$genome" = "$reference" ]
    then
        echo "Skipping reference genome itself: $genome"
        continue
    fi

    if ! is_fasta "$genome"
    then
        continue
    fi

    name=$(basename "$genome")
    name="${name%.*}"
    mkdir -p "${output_dir}/${name}"

    echo "Processing $genome ...."

    # Re-check index each genome (another run may have wiped a shared index before).
    ensure_reference_index "$ref_local"

    refbase=$(basename "$reference")
    refbase="${refbase%.*}"
    # Keep a short label from the reference name (last underscore token if present).
    if [[ "$refbase" == *_* ]]
    then
        refname="${refbase##*_}"
    else
        refname="$refbase"
    fi
    out_name="${name}_vs_${refname}"

    bwa mem -t "$cpus" -o "${output_dir}/${name}/${out_name}.sam" "$ref_local" "$genome"
    samtools view -Sb "${output_dir}/${name}/${out_name}.sam" > "${output_dir}/${name}/${out_name}.bam"
    rm -f "${output_dir}/${name}/${out_name}.sam"

    samtools sort "${output_dir}/${name}/${out_name}.bam" -o "${output_dir}/${name}/${out_name}_al.bam"
    rm -f "${output_dir}/${name}/${out_name}.bam"
    samtools index "${output_dir}/${name}/${out_name}_al.bam"

    bcftools mpileup -Ou -f "$ref_local" "${output_dir}/${name}/${out_name}_al.bam" --threads "$cpus" \
        | bcftools call -mv -Ov -o "${output_dir}/${name}/${out_name}.vcf"

    echo "Variant calls for $genome are ready at '${output_dir}/${name}'"
    processed=$((processed + 1))
done

if [ "$processed" -eq 0 ]
then
    echo "No query genomes left to compare after skipping the reference. Add other FASTA files or choose a different wild-type reference." >&2
    exit 1
fi

echo "SNP finder finished for ${processed} genome(s)."
