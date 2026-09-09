#!/usr/bin/env bash
set -euo pipefail

export BACTFLOW_IN_DOCKER=1
export BACTFLOW_MODULE="${BACTFLOW_MODULE:-assem}"
export BACTFLOW_HOST_PORT="${BACTFLOW_HOST_PORT:-5002}"
# Do not force BACTFLOW_NO_BROWSER — let Flask ask the host browser hook
# (or open locally). bactflow.sh may still set NO_BROWSER=1 itself.
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
