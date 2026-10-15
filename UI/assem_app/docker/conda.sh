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

# Real activate/deactivate in-shell (the /opt/conda/bin/conda stub cannot
# mutate the parent environment). Supports bactflow + bactflow_mlst sidecar.
conda() {
    local cmd="${1:-}"
    if [ "${cmd}" = "activate" ]; then
        local env_name="${2:-bactflow}"
        local prefix="/opt/conda/envs/${env_name}"
        if [ ! -d "${prefix}" ]; then
            echo "conda activate: unknown env '${env_name}'" >&2
            return 1
        fi
        # Drop previous env bin from PATH, then prepend the target.
        if [ -n "${CONDA_PREFIX:-}" ]; then
            PATH="$(printf '%s' "${PATH}" | awk -v p="${CONDA_PREFIX}/bin" -F: '{
                for (i=1;i<=NF;i++) if ($i != p && $i != "") out = (out ? out ":" : "") $i
                print out
            }')"
        fi
        export CONDA_PREFIX="${prefix}"
        export CONDA_DEFAULT_ENV="${env_name}"
        export CONDA_PYTHON_EXE="${prefix}/bin/python"
        export PATH="${prefix}/bin:/opt/conda/bin:${PATH}"
        if [ -d "${prefix}/lib" ]; then
            export LD_LIBRARY_PATH="${prefix}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
        fi
        # Sidecar mlst must not inherit the main env's Perl @INC.
        if [ "${env_name}" = "bactflow_mlst" ]; then
            unset PERL5LIB PERL_LOCAL_LIB_ROOT PERL_MB_OPT PERL_MM_OPT 2>/dev/null || true
        fi
        return 0
    fi
    if [ "${cmd}" = "deactivate" ]; then
        conda activate bactflow
        return $?
    fi
    /opt/conda/bin/conda "$@"
}
