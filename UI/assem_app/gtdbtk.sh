#!/bin/bash
# GTDB-Tk wrapper. Do not pip-upgrade gtdbtk: 2.7+ needs R232, 2.4.1–2.6.1 need R220/R226.

set -euo pipefail

genomes="./"
cpus=35
extension="fasta"
out_dir="./gtdbtk_out"
db_dir="${GTDBTK_DATA_PATH:-}"

display_help(){
    echo "Usage: $0 -g <genomes_directory> -c <cpus> -e <extension> -d <database directory> -o <output_directory>"
    echo "Options:"
    echo "  -g   Genomes directory to FASTA files"
    echo "  -c   Number of CPUs (default: 35)"
    echo "  -e   File extension without dot (default: fasta)"
    echo "  -d   Path to the GTDB-Tk reference data"
    echo "  -o   Output directory (default: gtdbtk_out)"
    exit 1
}

while getopts ":g:c:e:d:o:" opt
do
    case $opt in
        g) genomes="$OPTARG" ;;
        c) cpus="$OPTARG" ;;
        e) extension="$OPTARG" ;;
        d) db_dir="$OPTARG" ;;
        o) out_dir="$OPTARG" ;;
        \?) echo "Invalid option: -$OPTARG" >&2; exit 1 ;;
        :) echo "Option -$OPTARG requires an argument." >&2; exit 1 ;;
    esac
done

if [ "$#" -eq 0 ]
then
    display_help
fi

mkdir -p "${out_dir}"
extension="${extension#.}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bactflow
export PATH="${CONDA_PREFIX:-}/bin:${PATH}"
export GTDBTK_DATA_PATH="${db_dir}"

echo "Checking GTDB-Tk installations ..."

if ! command -v gtdbtk >/dev/null 2>&1; then
    echo "ERROR: gtdbtk is not on PATH." >&2
    echo "       conda activate bactflow && conda install -c bioconda 'gtdbtk=2.6.1'" >&2
    echo "       (2.6.1 matches GTDB R226. 2.7+ requires R232.)" >&2
    exit 1
fi

# pplacer ships a binary named guppy (tree conversion). This is NOT Oxford Nanopore Guppy.
if ! command -v pplacer >/dev/null 2>&1 || ! command -v guppy >/dev/null 2>&1; then
    echo "GTDB-Tk classify needs pplacer and its guppy helper (not ONT Guppy)."
    if command -v micromamba >/dev/null 2>&1; then
        micromamba install -y -p "${CONDA_PREFIX:-}" -c bioconda -c conda-forge pplacer || true
    elif command -v mamba >/dev/null 2>&1; then
        mamba install -y -c bioconda -c conda-forge pplacer || true
    fi
    hash -r 2>/dev/null || true
fi
if ! command -v pplacer >/dev/null 2>&1 || ! command -v guppy >/dev/null 2>&1; then
    echo "ERROR: pplacer/guppy is not on PATH." >&2
    echo "       Host: conda install -c bioconda pplacer" >&2
    echo "       Docker: rebuild the post-assembly image (keep bin/guppy)." >&2
    exit 1
fi
if ! command -v skani >/dev/null 2>&1; then
    echo "ERROR: skani is not on PATH (required by GTDB-Tk 2.4+ classify)." >&2
    echo "       conda install -c bioconda skani" >&2
    exit 1
fi
if ! command -v prodigal >/dev/null 2>&1 || ! command -v hmmsearch >/dev/null 2>&1; then
    echo "ERROR: prodigal and hmmer are required for gtdbtk identify." >&2
    exit 1
fi

echo "pplacer: $(command -v pplacer)"
echo "pplacer guppy: $(command -v guppy)"
echo "skani: $(command -v skani)"

tk_ver="$(gtdbtk -v 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)"
echo "I found gtdbtk version ${tk_ver:-unknown} in $(command -v gtdbtk)"

if [ -z "${db_dir}" ] || [ ! -d "${db_dir}" ]; then
    echo "ERROR: GTDB-Tk database directory is missing: ${db_dir:-<empty>}" >&2
    exit 1
fi
export GTDBTK_DATA_PATH="${db_dir}"
echo "GTDB-Tk database: ${db_dir}"

db_rel=""
if [ -f "${db_dir}/metadata/metadata.txt" ]; then
    db_rel="$(grep -Eo '[Rr]?[0-9]{3}' "${db_dir}/metadata/metadata.txt" | head -1 | tr -d 'Rr' || true)"
fi
if [ -z "${db_rel}" ]; then
    db_rel="$(basename "${db_dir}" | grep -Eo '[0-9]{3}' | head -1 || true)"
fi
if [ -n "${db_rel}" ]; then
    echo "Detected reference data release: R${db_rel}"
fi

# Official compatibility (https://ecogenomics.github.io/GTDBTk/installing/index.html):
#   R232: 2.7.0+
#   R226: 2.4.1 – 2.6.1
#   R220: 2.4.0 – 2.6.1
if [ -n "${tk_ver}" ] && [ -n "${db_rel}" ]; then
    tk_minor="${tk_ver#*.}"
    tk_minor="${tk_minor%%.*}"
    if [ "${tk_ver%%.*}" -ge 2 ] && [ "${tk_minor}" -ge 7 ] && [ "${db_rel}" != "232" ]; then
        echo "ERROR: GTDB-Tk ${tk_ver} requires reference data R232, but this database is R${db_rel}." >&2
        echo "       That is why classify started failing after the 2.7 upgrade." >&2
        echo "       Keep R226:  conda install -c bioconda 'gtdbtk=2.6.1'" >&2
        echo "       Or download R232 and point --gtdbtk_data_path at it." >&2
        exit 1
    fi
    if [ "${tk_ver%%.*}" -eq 2 ] && [ "${tk_minor}" -le 6 ] && [ "${db_rel}" = "232" ]; then
        echo "ERROR: GTDB-Tk ${tk_ver} cannot use R232. Install gtdbtk>=2.7.0 or use R226." >&2
        exit 1
    fi
fi

if [ -f "${genomes}" ]; then
    mkdir -p "${out_dir}/genomes_in"
    cp -f "${genomes}" "${out_dir}/genomes_in/"
    genomes="${out_dir}/genomes_in"
fi
if [ ! -d "${genomes}" ]; then
    echo "ERROR: genome directory not found: ${genomes}" >&2
    exit 1
fi

shopt -s nullglob
hits=("${genomes}"/*."${extension}")
if [ "${#hits[@]}" -eq 0 ]; then
    for cand in fasta fa fna; do
        hits=("${genomes}"/*."${cand}")
        if [ "${#hits[@]}" -gt 0 ]; then
            extension="${cand}"
            echo "Using detected genome extension: ${extension}"
            break
        fi
    done
fi
if [ "${#hits[@]}" -eq 0 ]; then
    echo "ERROR: No FASTA files in ${genomes}" >&2
    ls -la "${genomes}" >&2 || true
    exit 1
fi

echo "Executing gene calling..."
gtdbtk identify --genome_dir "${genomes}" --out_dir "${out_dir}/identify" --cpus "${cpus}" --extension "${extension}"

echo "Executing aligning..."
gtdbtk align --identify_dir "${out_dir}/identify" --out_dir "${out_dir}/align" --cpus "${cpus}"

echo "Executing classification..."
classify_help="$(gtdbtk classify --help 2>&1 || true)"
classify_extra=()
if echo "${classify_help}" | grep -q -- '--skip_ani_screen'; then
    classify_extra+=(--skip_ani_screen)
fi
# pplacer with many CPUs is OOM-killed: it maps ~80GB then forks, and Docker/cgroup
# kills the process ("Killed") during "Caching likelihood information on reference tree".
pplacer_cpus="${GTDBTK_PPLACER_CPUS:-1}"
if ! [[ "${pplacer_cpus}" =~ ^[0-9]+$ ]] || [ "${pplacer_cpus}" -lt 1 ]; then
    pplacer_cpus=1
fi
if echo "${classify_help}" | grep -q -- '--pplacer_cpus'; then
    classify_extra+=(--pplacer_cpus "${pplacer_cpus}")
fi
scratch_dir="${out_dir}/pplacer_scratch"
mkdir -p "${scratch_dir}"
if echo "${classify_help}" | grep -q -- '--scratch_dir'; then
    classify_extra+=(--scratch_dir "${scratch_dir}")
fi
echo "GTDB-Tk classify: --cpus ${cpus} --pplacer_cpus ${pplacer_cpus} --scratch_dir ${scratch_dir}"
echo "Using a pplacer scratch file so the ~80 GB reference-tree cache is not held in RAM."
gtdbtk classify --genome_dir "${genomes}" --align_dir "${out_dir}/align" --out_dir "${out_dir}/classify" -x "${extension}" --cpus "${cpus}" "${classify_extra[@]}"

shopt -s nullglob
summaries=("${out_dir}"/classify/*.summary.tsv)
if [ "${#summaries[@]}" -eq 0 ]; then
    echo "ERROR: GTDB-Tk classify finished without a summary TSV in ${out_dir}/classify" >&2
    ls -la "${out_dir}" "${out_dir}/classify" >&2 || true
    exit 1
fi
echo "GTDB-Tk summaries ready:"
ls -l "${summaries[@]}"
# Copy summaries next to the process (Nextflow output glob) and into GTDBTK_PUBLISH_DIR
# (the mounted results folder). Do not rely on publishDir of the huge identify/align tree.
for f in "${summaries[@]}"
do
    base="$(basename "$f")"
    case "${base}" in
        *markers*) continue ;;
    esac
    cp -f "$f" "${out_dir}/classify/${base}"
    cp -f "$f" "./${base}"
    if [ -n "${GTDBTK_PUBLISH_DIR:-}" ]
    then
        mkdir -p "${GTDBTK_PUBLISH_DIR}"
        cp -f "$f" "${GTDBTK_PUBLISH_DIR}/${base}"
    fi
done
if [ -n "${GTDBTK_PUBLISH_DIR:-}" ]
then
    echo "Published GTDB-Tk summaries to ${GTDBTK_PUBLISH_DIR}"
    ls -l "${GTDBTK_PUBLISH_DIR}"/*.summary.tsv
    shopt -s nullglob
    for tree in "${out_dir}"/classify/*.classify.tree "${out_dir}"/classify/*.decorated.tree
    do
        [ -f "${tree}" ] || continue
        cp -f "${tree}" "${GTDBTK_PUBLISH_DIR}/$(basename "${tree}")"
    done
fi
echo "GTDB-Tk output directory: $(pwd)/${out_dir}"
ls -l ./*.summary.tsv
