#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=bactflow_runtime.sh
source "${SCRIPT_DIR}/bactflow_runtime.sh"

bactflow_setup_nextflow_env() {
    if command -v conda >/dev/null 2>&1; then
        # shellcheck disable=SC1091
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate bactflow
    fi

    bactflow_prepare_nextflow
}
