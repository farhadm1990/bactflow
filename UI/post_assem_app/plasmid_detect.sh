#!/bin/bash
# geNomad plasmid + virus (MGE) detection on assembled FASTA files.
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

mkdir -p "${out_dir}/plasmids" "${out_dir}/viruses" "${out_dir}/summaries" "${out_dir}/work"
: > "${out_dir}/genomes_scanned.txt"
echo "geNomad $(genomad --version 2>/dev/null | head -1 || echo unknown)"
echo "Database: ${db_dir}"
echo "Scanning ${#fastas[@]} genome FASTA(s) for plasmids and viruses..."

gn_extra=()
if genomad end-to-end --help 2>&1 | grep -q -- '--cleanup'
then
    gn_extra+=(--cleanup)
fi

for fa in "${fastas[@]}"
do
    stem="$(basename "${fa}")"
    stem="${stem%.*}"
    echo "${stem}" >> "${out_dir}/genomes_scanned.txt"
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

    psum="$(find "${work}" -type f -name '*_plasmid_summary.tsv' | head -n 1 || true)"
    if [ -n "${psum}" ] && [ -f "${psum}" ]
    then
        cp -f "${psum}" "${out_dir}/summaries/${stem}.plasmid_summary.tsv"
    fi
    pfna="$(find "${work}" -type f -name '*_plasmid.fna' | head -n 1 || true)"
    if [ -n "${pfna}" ] && [ -s "${pfna}" ]
    then
        cp -f "${pfna}" "${out_dir}/plasmids/${stem}_plasmids.fna"
    fi

    vsum="$(find "${work}" -type f -name '*_virus_summary.tsv' | head -n 1 || true)"
    if [ -n "${vsum}" ] && [ -f "${vsum}" ]
    then
        cp -f "${vsum}" "${out_dir}/summaries/${stem}.virus_summary.tsv"
    fi
    vfna="$(find "${work}" -type f -name '*_virus.fna' | head -n 1 || true)"
    if [ -n "${vfna}" ] && [ -s "${vfna}" ]
    then
        cp -f "${vfna}" "${out_dir}/viruses/${stem}_viruses.fna"
    fi
done

python3 - "${out_dir}" <<'PY'
import csv
import os
import sys

out_dir = sys.argv[1]
sum_dir = os.path.join(out_dir, "summaries")


def merge_summaries(suffix, dest_name, default_header):
    dest = os.path.join(out_dir, dest_name)
    rows = []
    header = []
    if os.path.isdir(sum_dir):
        for name in sorted(os.listdir(sum_dir)):
            if not name.endswith(suffix):
                continue
            path = os.path.join(sum_dir, name)
            genome = name[: -len(suffix)]
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
                    rows.append({"Genome": genome, **rec})
    if not header:
        header = default_header
    with open(dest, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} row(s) to {dest}")
    return dest


merge_summaries(
    ".plasmid_summary.tsv",
    "plasmid_summary.tsv",
    ["Genome", "seq_name", "length", "topology", "plasmid_score", "n_hallmarks", "conjugation_genes", "amr_genes"],
)
merge_summaries(
    ".virus_summary.tsv",
    "virus_summary.tsv",
    ["Genome", "seq_name", "length", "topology", "coordinates", "virus_score", "n_hallmarks", "marker_enrichment", "taxonomy"],
)
PY

echo "geNomad finished."
echo "  Plasmids: ${out_dir}/plasmid_summary.tsv  FASTA: ${out_dir}/plasmids"
echo "  Viruses:  ${out_dir}/virus_summary.tsv   FASTA: ${out_dir}/viruses"
ls -l "${out_dir}/plasmid_summary.tsv" "${out_dir}/virus_summary.tsv" "${out_dir}/plasmids" "${out_dir}/viruses" 2>/dev/null || true
