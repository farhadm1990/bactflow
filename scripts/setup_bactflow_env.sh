#!/usr/bin/env bash
# Create or update the bactflow Conda env from config.yml.
# Does not download GTDB-Tk / CheckM / Bakta / geNomad databases — those are
# supplied later in the UI.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_FILE="${1:-${ROOT_DIR}/config.yml}"
ENV_NAME="${2:-bactflow}"

cd "${ROOT_DIR}"

if ! command -v mamba >/dev/null 2>&1 && ! command -v conda >/dev/null 2>&1; then
    echo "conda/mamba not found in PATH" >&2
    exit 1
fi

# CheckM's conda post-link.sh wget's a DB unless this marker already exists.
# Leave a placeholder so env create stays offline; the UI still asks for --checkm_db.
export CHECKM_DATA_DIR="${CHECKM_DATA_DIR:-${HOME}/.bactflow/checkm_skip_postlink}"
mkdir -p "${CHECKM_DATA_DIR}/genome_tree"
if [[ ! -e "${CHECKM_DATA_DIR}/genome_tree/genome_tree.derep.txt" ]]; then
    : > "${CHECKM_DATA_DIR}/genome_tree/genome_tree.derep.txt"
fi

# Do not let AMRFinder / GTDB-Tk post-install hooks fetch databases.
export AMRFINDER_IGNORE_MISSING_DB=1
unset GTDBTK_DATA_PATH || true

echo "Creating/updating Conda env '${ENV_NAME}' from ${ENV_FILE}"
echo "Databases are not downloaded here. Point Bakta / GTDB-Tk / CheckM / geNomad at local DBs in the UI."

# conda activate/deactivate scripts are not `set -u` safe (MKL backup vars).
set +u
if command -v mamba >/dev/null 2>&1; then
    mamba env create -f "${ENV_FILE}" -y || mamba env update -n "${ENV_NAME}" -f "${ENV_FILE}" -y
else
    conda env create -f "${ENV_FILE}" -y || conda env update -n "${ENV_NAME}" -f "${ENV_FILE}" -y
fi
set -u

bash "${ROOT_DIR}/scripts/install_bactflow_hooks.sh" "${ENV_NAME}"

# PubMLST (mlst) cannot share the main env (perl-bioperl vs blast/samtools).
MLST_ENV_FILE="${ROOT_DIR}/UI/bactflow_app/environment.mlst.yml"
MLST_ENV_NAME="bactflow_mlst"
if [[ -f "${MLST_ENV_FILE}" ]]; then
    echo "Creating/updating sidecar Conda env '${MLST_ENV_NAME}' from ${MLST_ENV_FILE}"
    set +u
    if command -v mamba >/dev/null 2>&1; then
        mamba env create -f "${MLST_ENV_FILE}" -y || mamba env update -n "${MLST_ENV_NAME}" -f "${MLST_ENV_FILE}" -y
    else
        conda env create -f "${MLST_ENV_FILE}" -y || conda env update -n "${MLST_ENV_NAME}" -f "${MLST_ENV_FILE}" -y
    fi
    set -u
fi

echo ""
echo "Done. Use: conda activate ${ENV_NAME}"
echo "MLST runs from sidecar env '${MLST_ENV_NAME}' (see mlst_run.sh)."
echo "No reference databases were downloaded. Set them in the UI when you run a module that needs them."
