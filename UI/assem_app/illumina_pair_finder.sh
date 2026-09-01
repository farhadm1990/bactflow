#!/usr/bin/env bash
set -euo pipefail

fastq_dir=""
extension=".fastq.gz"

usage() {
    echo "Usage: $0 -g <fastq_directory> -e <extension>"
    exit 1
}

while getopts "g:e:" opt; do
    case "${opt}" in
        g) fastq_dir="${OPTARG}" ;;
        e) extension="${OPTARG}" ;;
        *) usage ;;
    esac
done

if [[ -z "${fastq_dir}" ]]; then
    usage
fi

fastq_dir=$(realpath "${fastq_dir}")
manifest=$(mktemp)

printf "sample\tr1\tr2\n" > "${manifest}"

stage_pair() {
    local sample="$1"
    local r1_in="$2"
    local r2_in="$3"
    printf "%s\t%s\t%s\n" "${sample}" "${r1_in}" "${r2_in}" >> "${manifest}"
}

concat_pair_in_dir() {
    local sample_dir="$1"
    local sample
    sample=$(basename "${sample_dir}")
    local r1_files=()
    local r2_files=()

    if ls "${sample_dir}"/*_R1*"${extension}" 1> /dev/null 2>&1; then
        r1_files=("${sample_dir}"/*_R1*"${extension}")
        r2_files=("${sample_dir}"/*_R2*"${extension}")
    elif ls "${sample_dir}"/*_1"${extension}" 1> /dev/null 2>&1; then
        r1_files=("${sample_dir}"/*_1"${extension}")
        r2_files=("${sample_dir}"/*_2"${extension}")
    else
        return 1
    fi

    stage_pair "${sample}" "${r1_files[*]}" "${r2_files[*]}"
}

for sample_dir in "${fastq_dir}"/*; do
    if [[ -d "${sample_dir}" ]]; then
        sample_name=$(basename "${sample_dir}")
        if [[ "${sample_name}" != "pooled" && "${sample_name}" != "paired" && "${sample_name}" != "pe_pairs" ]]; then
            concat_pair_in_dir "${sample_dir}" || true
        fi
    fi
done

for r1_file in "${fastq_dir}"/*_R1*"${extension}" "${fastq_dir}"/*_1"${extension}"; do
    if [[ ! -f "${r1_file}" ]]; then
        continue
    fi

    sample=$(basename "${r1_file}" | sed -E 's/_R1.*$//; s/_1\.fastq(\.gz)?$//')
    if [[ "${extension}" == ".fastq.gz" ]]; then
        r2_file=$(echo "${r1_file}" | sed -E 's/_R1/_R2/; s/_1\.fastq\.gz/_2.fastq.gz/')
    else
        r2_file=$(echo "${r1_file}" | sed -E 's/_R1/_R2/; s/_1\.fastq/_2.fastq/')
    fi

    if [[ -f "${r2_file}" ]]; then
        stage_pair "${sample}" "${r1_file}" "${r2_file}"
    fi
done

if [[ "$(wc -l < "${manifest}")" -le 1 ]]; then
    echo "No Illumina paired-end reads found in ${fastq_dir}" >&2
    rm -f "${manifest}"
    exit 1
fi

cat "${manifest}"
rm -f "${manifest}"
