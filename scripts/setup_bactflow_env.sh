#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_FILE="${1:-${ROOT_DIR}/config.yml}"
ENV_NAME="${2:-bactflow}"

cd "${ROOT_DIR}"

if ! command -v mamba >/dev/null 2>&1 && ! command -v conda >/dev/null 2>&1; then
    echo "conda/mamba not found in PATH" >&2
    exit 1
fi

if command -v mamba >/dev/null 2>&1; then
    mamba env create -f "${ENV_FILE}" -y || mamba env update -n "${ENV_NAME}" -f "${ENV_FILE}" -y
else
    conda env create -f "${ENV_FILE}" -y || conda env update -n "${ENV_NAME}" -f "${ENV_FILE}" -y
fi

bash "${ROOT_DIR}/scripts/install_bactflow_hooks.sh" "${ENV_NAME}"

echo ""
echo "Done. Use: conda activate ${ENV_NAME}"
