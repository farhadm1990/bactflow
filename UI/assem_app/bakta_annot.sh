#!/bin/bash

genomes="./"
cpus=35
out_dir="bakta_out"
db_dir=""

display_help(){
    echo "Usage: $0 -g <genomes_directory> -c <cpus> -o <output_directory>" #$0 is the name of the script
    echo "Options:"
    echo "  -g   Genomes directory to fasta files (default: ./)"
    echo "  -c   Number of CPUs (default: 35)"
    echo "  -d   Directory to bakta database"
    echo "  -o   Output directory (default: out_dir)"
    exit 0
}

while getopts ":g:c:d:o:" opt
do 
    case $opt in
        g) 
            genomes="$OPTARG"
            ;;
        c)
            cpus="$OPTARG"
            ;;
        d)
            db_dir="$OPTARG"
            ;;
    
        o) 
            out_dir="$OPTARG"
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

pick_bakta_bin() {
    if [ -n "${CONDA_PREFIX:-}" ] && [ -x "${CONDA_PREFIX}/bin/bakta" ]; then
        echo "${CONDA_PREFIX}/bin/bakta"
        return 0
    fi
    if command -v bakta >/dev/null 2>&1; then
        command -v bakta
        return 0
    fi
    return 1
}

BAKTA_BIN="$(pick_bakta_bin)" || {
    echo "ERROR: bakta is not installed." >&2
    echo "       conda activate bactflow && conda install -c bioconda -c conda-forge 'bakta>=1.11'" >&2
    exit 1
}

if [[ "$BAKTA_BIN" == *"/.local/"* ]]; then
    echo "WARNING: Using user pip install ($BAKTA_BIN). Prefer conda: conda install -c bioconda -c conda-forge 'bakta>=1.11'" >&2
fi

if [ -z "$db_dir" ] || [ ! -d "$db_dir" ]; then
    echo "ERROR: Bakta database directory not found: ${db_dir:-<empty>}" >&2
    exit 1
fi

if ! command -v "$BAKTA_BIN" >/dev/null 2>&1 && [ ! -x "$BAKTA_BIN" ]; then
    echo "ERROR: bakta binary not executable: $BAKTA_BIN" >&2
    exit 1
fi

bakta_ver="$("$BAKTA_BIN" --version 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)"
bakta_major="${bakta_ver%%.*}"
bakta_minor="$(echo "$bakta_ver" | cut -d. -f2)"

db_version_file=""
for candidate in "${db_dir}/db.version" "${db_dir}/version.json"; do
    if [ -f "$candidate" ]; then
        db_version_file="$candidate"
        break
    fi
done

if [ -n "$db_version_file" ]; then
    if [[ "$db_version_file" == *.json ]]; then
        db_major="$(python3 -c "import json; d=json.load(open('${db_version_file}')); print(str(d.get('major', d.get('version', '0'))).split('.')[0])" 2>/dev/null || echo "0")"
    else
        db_major="$(cut -d. -f1 < "$db_version_file")"
    fi

    if [ "$db_major" -ge 6 ] 2>/dev/null; then
        if [ "$bakta_major" -lt 1 ] || { [ "$bakta_major" -eq 1 ] && [ "$bakta_minor" -lt 11 ]; }; then
            echo "ERROR: Bakta ${bakta_ver:-unknown} requires database v5.x, but v${db_major}.x was detected at:" >&2
            echo "       ${db_dir}" >&2
            echo "       Update Bakta: conda install -c bioconda -c conda-forge 'bakta>=1.11'" >&2
            echo "       Or download a v5.x database from https://doi.org/10.5281/zenodo.4247253" >&2
            exit 1
        fi
    fi
fi

if [ -n "${CONDA_PREFIX:-}" ]; then
    export PATH="${CONDA_PREFIX}/bin:${PATH}"
fi

require_bakta_tool() {
    local tool="$1"
    local hint="$2"
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: ${tool} not found. ${hint}" >&2
        echo "       bash scripts/install_bakta_deps.sh bactflow" >&2
        exit 1
    fi
}

require_bakta_tool "tRNAscan-SE" "Bakta requires tRNAscan-SE v2.0.12."
if ! tRNAscan-SE -h 2>&1 | grep -qE 'tRNAscan-SE 2\.0\.(1[12]|[0-9]{2})'; then
    echo "ERROR: Wrong tRNAscan-SE version (Bakta needs v2.0.12)." >&2
    echo "       Found: $(tRNAscan-SE -h 2>&1 | head -1)" >&2
    echo "       Fix: bash scripts/install_bakta_deps.sh bactflow" >&2
    exit 1
fi
require_bakta_tool "cmscan" "Bakta requires Infernal cmscan v1.1.5."
if ! cmscan -h 2>&1 | grep -qE 'INFERNAL 1\.1\.5'; then
    echo "ERROR: Wrong CMscan/Infernal version (Bakta needs v1.1.5)." >&2
    echo "       Found: $(cmscan -h 2>&1 | sed -n '2p')" >&2
    echo "       Fix: bash scripts/install_bakta_deps.sh bactflow" >&2
    exit 1
fi
require_bakta_tool "aragorn" "Required by Bakta for tmRNA annotation."

cpus=$(($cpus))
shopt -s nullglob
fasta_files=("$genomes"/*.fasta)
if [ "${#fasta_files[@]}" -eq 0 ]; then
    echo "ERROR: No .fasta files found in ${genomes}" >&2
    exit 1
fi

for file in "${fasta_files[@]}"
do 
    striped="$(basename "$file" | cut -d'.' -f1)"
    dest="${out_dir}/${striped}_bakta"

    if ! "$BAKTA_BIN" --db "${db_dir}" "${file}" --prefix "$striped" --output "$dest" --force --threads "${cpus}"; then
        echo "ERROR: Bakta failed for ${file}" >&2
        exit 1
    fi
done
