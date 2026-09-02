#!/usr/bin/env bash
# Install pinned Bakta binary dependencies (conda solve fails when prokka shares the env).
set -euo pipefail

ENV_NAME="${1:-bactflow}"
BIOCONDA="https://conda.anaconda.org/bioconda/linux-64"

install_pkg() {
    local pkg="$1"
    local url="${BIOCONDA}/${pkg}"
    local dest="${tmpdir}/${pkg}"
    echo "Installing ${pkg} into ${CONDA_PREFIX}..."
    curl -fsSL -o "${dest}" "${url}"
    conda install -y --no-deps -p "${CONDA_PREFIX}" "${dest}"
}

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

if ! command -v tRNAscan-SE >/dev/null 2>&1 || ! tRNAscan-SE -h 2>&1 | grep -q "tRNAscan-SE 2.0.12"; then
    install_pkg "trnascan-se-2.0.12-pl5321h7b50bb2_2.tar.bz2"
fi

if ! command -v cmscan >/dev/null 2>&1 || ! cmscan -h 2>&1 | grep -q "INFERNAL 1.1.5"; then
    install_pkg "infernal-1.1.5-pl5321h7b50bb2_4.tar.bz2"
fi

if ! command -v aragorn >/dev/null 2>&1; then
    echo "Installing aragorn via conda (with deps)..."
    conda install -y -c bioconda aragorn || conda install -y -c bioconda --no-deps aragorn
fi

echo "Bakta dependency check:"
echo "  $(tRNAscan-SE -h 2>&1 | head -1)"
echo "  $(cmscan -h 2>&1 | sed -n '2p')"
command -v aragorn >/dev/null && echo "  aragorn: $(aragorn -h 2>&1 | head -1 || true)"
