#!/usr/bin/env bash
set -euo pipefail

export BACTFLOW_IN_DOCKER=1
export BACTFLOW_NO_BROWSER=1
export PYTHONUNBUFFERED=1
export CONDA_PREFIX="${CONDA_PREFIX:-/opt/conda/envs/bactflow}"
export CONDA_DEFAULT_ENV="${CONDA_DEFAULT_ENV:-bactflow}"
export CONDA_SHLVL=1
export PATH="${CONDA_PREFIX}/bin:/opt/conda/bin:/usr/local/bin:${PATH}"
if [ -d "${CONDA_PREFIX}/lib" ]; then
    export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi

if [ -d /usr/lib/jvm/bactflow-java ]; then
    export JAVA_HOME=/usr/lib/jvm/bactflow-java
    export BACTFLOW_JAVA_HOME="${JAVA_HOME}"
    export PATH="${JAVA_HOME}/bin:${PATH}"
fi

if [ -f /etc/profile.d/bactflow.sh ]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/bactflow.sh
elif [ -f /opt/conda/etc/profile.d/conda.sh ]; then
    # shellcheck disable=SC1091
    source /opt/conda/etc/profile.d/conda.sh
fi

if [ "$#" -eq 0 ]; then
    echo "No command given to Docker entrypoint" >&2
    exit 1
fi

exec "$@"
