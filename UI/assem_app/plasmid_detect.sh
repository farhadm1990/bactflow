#!/bin/bash
# geNomad plasmid detection on assembled FASTA files.
# Database: genomad download-database /path/to/genomad_db

set -euo pipefail

genomes="./"
cpus=8
db_dir=""
out_dir="plasmid_out"

display_help(){
    echo "Usage: $0 -g <genomes_directory> -c <cpus> -d <genomad_db> -o <output_directory>"
    echo "Options:"
    echo "  -g   Directory of FASTA genomes (default: ./)"
    echo "  -c   CPUs (default: 8)"
    echo "  -d   geNomad database directory (from: genomad download-database DIR)"
    echo "  -o   Output directory (default: plasmid_out)"
    exit 1
}

while getopts ":g:c:d:o:" opt
do
    case $opt in
        g) genomes="$OPTARG" ;;
        c) cpus="$OPTARG" ;;
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

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate bactflow
export PATH="${CONDA_PREFIX:-}/bin:${PATH}"

if ! command -v genomad >/dev/null 2>&1
then
    echo "ERROR: genomad is not on PATH." >&2
    echo "       conda install -c bioconda genomad" >&2
    echo "       Then: genomad download-database /path/to/genomad_db" >&2
    exit 1
fi

if [ -z "${db_dir}" ]
then
    echo "ERROR: geNomad database directory is missing." >&2
    echo "       genomad download-database /path/to/genomad_db" >&2
    exit 1
fi
if [ -d "${db_dir}/genomad_db" ]
then
    db_dir="${db_dir}/genomad_db"
fi
if [ ! -d "${db_dir}" ]
then
    echo "ERROR: geNomad database directory not found: ${db_dir}" >&2
    exit 1
fi

if [ -f "${genomes}" ]
then
    mkdir -p genomes_in
    cp -f "${genomes}" genomes_in/
    genomes="genomes_in"
fi
if [ ! -d "${genomes}" ]
then
    echo "ERROR: genome directory not found: ${genomes}" >&2
    exit 1
fi

shopt -s nullglob
fastas=("${genomes}"/*.fasta "${genomes}"/*.fa "${genomes}"/*.fna)
if [ "${#fastas[@]}" -eq 0 ]
then
    echo "ERROR: No FASTA files in ${genomes}" >&2
    ls -la "${genomes}" >&2 || true
    exit 1
fi

mkdir -p "${out_dir}/plasmids" "${out_dir}/summaries" "${out_dir}/work"
echo "geNomad $(genomad --version 2>/dev/null | head -1 || echo unknown)"
echo "Database: ${db_dir}"
echo "Scanning ${#fastas[@]} genome FASTA(s) for plasmids..."

gn_extra=()
if genomad end-to-end --help 2>&1 | grep -q -- '--cleanup'
then
    gn_extra+=(--cleanup)
fi

for fa in "${fastas[@]}"
do
    stem="$(basename "${fa}")"
    stem="${stem%.*}"
    work="${out_dir}/work/${stem}"
    mkdir -p "${work}"
    extra=("${gn_extra[@]}")
    size="$(wc -c < "${fa}" | tr -d ' ')"
    if [ "${size}" -gt 80000000 ] && genomad end-to-end --help 2>&1 | grep -q -- '--splits'
    then
        extra+=(--splits 8)
    fi
    echo "geNomad: ${stem}"
    genomad end-to-end "${fa}" "${work}" "${db_dir}" --threads "${cpus}" "${extra[@]}"
    sum="$(find "${work}" -type f -name '*_plasmid_summary.tsv' | head -n 1 || true)"
    if [ -n "${sum}" ] && [ -f "${sum}" ]
    then
        cp -f "${sum}" "${out_dir}/summaries/${stem}.plasmid_summary.tsv"
    fi
    fna="$(find "${work}" -type f -name '*_plasmid.fna' | head -n 1 || true)"
    if [ -n "${fna}" ] && [ -s "${fna}" ]
    then
        cp -f "${fna}" "${out_dir}/plasmids/${stem}_plasmids.fna"
    fi
done

python3 - "${out_dir}" <<'PY'
import csv
import os
import sys

out_dir = sys.argv[1]
sum_dir = os.path.join(out_dir, "summaries")
dest = os.path.join(out_dir, "plasmid_summary.tsv")
rows = []
header = []
if os.path.isdir(sum_dir):
    for name in sorted(os.listdir(sum_dir)):
        if not name.endswith(".plasmid_summary.tsv"):
            continue
        path = os.path.join(sum_dir, name)
        genome = name[: -len(".plasmid_summary.tsv")]
        with open(path, encoding="utf-8", errors="replace") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if not reader.fieldnames:
                continue
            if not header:
                header = ["Genome"] + list(reader.fieldnames)
            for rec in reader:
                if not rec:
                    continue
                contig = (rec.get("seq_name") or rec.get("Contig") or "").strip()
                if not contig:
                    continue
                rec = {"Genome": genome, **rec}
                rows.append(rec)
if not header:
    header = ["Genome", "seq_name", "length", "topology", "plasmid_score", "n_hallmarks", "conjugation_genes", "amr_genes"]
with open(dest, "w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    writer.writerows(rows)
print(f"Wrote {len(rows)} plasmid contig(s) to {dest}")
PY

echo "Plasmid detection finished. Summary: ${out_dir}/plasmid_summary.tsv"
ls -l "${out_dir}/plasmid_summary.tsv" "${out_dir}/plasmids" 2>/dev/null || true
