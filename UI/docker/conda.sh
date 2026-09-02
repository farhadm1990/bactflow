# Conda-compatible activate script for packed Docker images.
# Existing BactFlow scripts call: source "$(conda info --base)/etc/profile.d/conda.sh"
export CONDA_PREFIX="${CONDA_PREFIX:-/opt/conda/envs/bactflow}"
export CONDA_DEFAULT_ENV="${CONDA_DEFAULT_ENV:-bactflow}"
export CONDA_SHLVL="${CONDA_SHLVL:-1}"
export CONDA_EXE="${CONDA_EXE:-/opt/conda/bin/conda}"
export CONDA_PYTHON_EXE="${CONDA_PREFIX}/bin/python"
export PATH="${CONDA_PREFIX}/bin:/opt/conda/bin:${PATH}"
if [ -d "${CONDA_PREFIX}/lib" ]; then
    export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
fi

conda() {
    /opt/conda/bin/conda "$@"
}
