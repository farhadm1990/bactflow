#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${1:-bactflow}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNTIME_SH="${SCRIPT_DIR}/bactflow_runtime.sh"
HOOK_NAME="zz-bactflow-runtime.sh"

if ! command -v conda >/dev/null 2>&1; then
    echo "conda not found in PATH" >&2
    exit 1
fi

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"

if ! conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
    echo "Conda env '${ENV_NAME}' not found. Create it first: bash scripts/setup_bactflow_env.sh" >&2
    exit 1
fi

CONDA_PREFIX="$(conda info --envs | awk -v env="${ENV_NAME}" '$1 == env {print $NF; exit}')"
if [ -z "${CONDA_PREFIX}" ] || [ ! -d "${CONDA_PREFIX}" ]; then
    echo "Could not resolve CONDA_PREFIX for env '${ENV_NAME}'" >&2
    exit 1
fi

HOOK_DIR="${CONDA_PREFIX}/etc/conda/activate.d"
mkdir -p "${HOOK_DIR}"
rm -f "${HOOK_DIR}/bactflow-nextflow-java.sh" "${HOOK_DIR}/bactflow-runtime.sh"
cp "${RUNTIME_SH}" "${HOOK_DIR}/${HOOK_NAME}"
chmod +x "${HOOK_DIR}/${HOOK_NAME}"

echo "Installed BactFlow runtime hook: ${HOOK_DIR}/${HOOK_NAME}"
echo "Re-activate with: conda deactivate && conda activate ${ENV_NAME}"

conda activate "${ENV_NAME}"
if bactflow_prepare_nextflow; then
    echo "OK: Java=$JAVA_CMD"
    echo "OK: Nextflow=$BACTFLOW_NEXTFLOW_BIN"
    "${BACTFLOW_NEXTFLOW_BIN}" -version 2>&1 | head -n1
else
    echo "Hook installed but verification failed — check messages above." >&2
    exit 1
fi
