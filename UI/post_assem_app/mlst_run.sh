#!/usr/bin/env bash
# Multi-locus sequence typing (Torsten Seemann mlst / PubMLST) on assembled FASTAs.
set -euo pipefail

genomes="./"
cpus=8
out_dir="mlst_out"

display_help(){
    echo "Usage: $0 -g <genomes_directory> -c <cpus> -o <output_directory>"
    echo "Options:"
    echo "  -g   Directory of FASTA genomes (default: ./)"
    echo "  -c   CPUs (default: 8)"
    echo "  -o   Output directory (default: mlst_out)"
    exit 1
}

while getopts ":g:c:o:" opt
do
    case $opt in
        g) genomes="$OPTARG" ;;
        c) cpus="$OPTARG" ;;
        o) out_dir="$OPTARG" ;;
        \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
    esac
done

if [ "$#" -eq 0 ]
then
    display_help
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
_conda_base="$(conda info --base)"
# Prefer bactflow_mlst sidecar (Docker + host): mlst's perl-bioperl stack
# cannot share an env with infernal/samtools/spades pins in bactflow.
if [ -x "${_conda_base}/envs/bactflow_mlst/bin/mlst" ]
then
    conda activate bactflow_mlst
elif [ -x "${_conda_base}/envs/bactflow/bin/mlst" ]
then
    conda activate bactflow
elif conda env list | awk '{print $1}' | grep -qx 'bactflow_mlst'
then
    conda activate bactflow_mlst
elif conda env list | awk '{print $1}' | grep -qx 'bactflow'
then
    conda activate bactflow
else
    echo "ERROR: neither bactflow_mlst nor bactflow conda env found." >&2
    exit 1
fi
export PATH="${CONDA_PREFIX:-}/bin:${PATH}"

if ! command -v mlst >/dev/null 2>&1
then
    echo "ERROR: mlst is not on PATH. Install the bioconda package 'mlst'." >&2
    echo "       conda create -n bactflow_mlst -c bioconda -c conda-forge mlst=2.23.0" >&2
    echo "       (or: conda install -n bactflow -c bioconda -c conda-forge mlst=2.23.0)" >&2
    exit 1
fi

if [ ! -d "${genomes}" ]
then
    echo "ERROR: genomes directory missing: ${genomes}" >&2
    exit 1
fi

# Nextflow often stages input dirs as symlinks. GNU find without -L / a
# trailing slash will not descend into a symlink-to-directory start path.
shopt -s nullglob
genomes_resolved="$(readlink -f "${genomes}" 2>/dev/null || echo "${genomes}")"
# Unique by basename stem so symlink + realpath double hits collapse to one genome.
mapfile -t fasta_files < <(
    {
        find -L "${genomes}" -maxdepth 1 -type f \( \
            -iname '*.fasta' -o -iname '*.fa' -o -iname '*.fna' \
        \) 2>/dev/null
        printf '%s\n' \
            "${genomes}"/*.fasta "${genomes}"/*.fa "${genomes}"/*.fna \
            "${genomes}"/*.FASTA "${genomes}"/*.FA "${genomes}"/*.FNA
    } | while IFS= read -r f
    do
        [ -n "$f" ] && [ -f "$f" ] || continue
        real="$(readlink -f "$f" 2>/dev/null || echo "$f")"
        base="$(basename "$real")"
        stem="${base%.*}"
        printf '%s\t%s\n' "${stem,,}" "$real"
    done | awk -F'\t' 'NF==2 && !seen[$1]++ {print $2}' | sort
)

if [ "${#fasta_files[@]}" -eq 0 ]
then
    echo "ERROR: No FASTA files in ${genomes}" >&2
    echo "       resolved: ${genomes_resolved}" >&2
    ls -la "${genomes}" "${genomes_resolved}" 2>&1 | head -40 >&2 || true
    exit 1
fi

mkdir -p "${out_dir}"
echo "mlst $(mlst --version 2>&1 | head -1 || echo unknown)"
echo "Typing ${#fasta_files[@]} genome(s) with ${cpus} thread(s)..."

# mlst writes NO header by default: FILE SCHEME ST allele...
# Capture raw rows, then prepend a header so the UI/report can parse it.
tmp_raw="${out_dir}/.mlst.raw.tsv"
mlst --threads "${cpus}" --nopath "${fasta_files[@]}" > "${tmp_raw}"

if [ ! -s "${tmp_raw}" ]
then
    echo "ERROR: mlst produced an empty table." >&2
    exit 1
fi

# Build header from the widest allele column count among rows.
max_alleles="$(awk -F'\t' '{n=NF-3; if(n>m) m=n} END{print m+0}' "${tmp_raw}")"
{
    printf 'FILE\tSCHEME\tST'
    if [ "${max_alleles}" -gt 0 ]
    then
        i=1
        while [ "${i}" -le "${max_alleles}" ]
        do
            printf '\tallele_%d' "${i}"
            i=$((i + 1))
        done
    fi
    printf '\n'
    # Prefer gene names from the first fully typed row (scheme != -).
    # Fall back to raw columns when schemes differ across samples.
    cat "${tmp_raw}"
} > "${out_dir}/mlst.tsv"
rm -f "${tmp_raw}"

# Compact summary: genome, scheme, ST (skip our header row).
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "genome","scheme","ST"; next} {print $1,$2,$3}' \
    "${out_dir}/mlst.tsv" > "${out_dir}/mlst_summary.tsv"

printf '%s\n' "${fasta_files[@]}" | xargs -n1 basename > "${out_dir}/genomes_scanned.txt"
echo "MLST written to ${out_dir}/mlst.tsv"
wc -l "${out_dir}/mlst.tsv"
