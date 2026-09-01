#!/usr/bin/env bash

# BactFlow runtime: compatible Java (17-24) + Nextflow (<26, pinned 24.10.2).
# Sourced by conda activate hooks, Docker images, and Flask UI subprocess wrappers.

BACTFLOW_NXF_PIN="${BACTFLOW_NXF_VER:-24.10.2}"
BACTFLOW_NXF_MAX_MAJOR=25

bactflow_java_major() {
    local ver major

    ver="$("$1" -version 2>&1 | awk -F '"' '/version/ {print $2; exit}')"
    major="$(echo "${ver}" | grep -oE '^[0-9]+' | head -n1)"

    if [ "${major}" = "1" ]; then
        major="$(echo "${ver}" | cut -d. -f2)"
    fi

    echo "${major}"
}

bactflow_java_is_nextflow_compatible() {
    local major

    major="$(bactflow_java_major "$1")"
    [ -n "${major}" ] && [ "${major}" -ge 17 ] && [ "${major}" -le 24 ]
}

bactflow_java_in_conda_env() {
    local resolved

    resolved="$(readlink -f "$1" 2>/dev/null || echo "$1")"
    case "${resolved}" in
        */envs/*|*/miniconda3/envs/*|*/anaconda3/envs/*) return 0 ;;
    esac
    return 1
}

bactflow_use_java_bin() {
    local java_bin="$1" home

    [ -n "${java_bin}" ] && [ -x "${java_bin}" ] || return 1
    bactflow_java_is_nextflow_compatible "${java_bin}" || return 1

    java_bin="$(readlink -f "${java_bin}" 2>/dev/null || echo "${java_bin}")"
    case "${java_bin}" in
        */lib/jvm/bin/java) home="${java_bin%/bin/java}" ;;
        *) home="$(dirname "$(dirname "${java_bin}")")" ;;
    esac
    export JAVA_HOME="${home}"
    export PATH="${JAVA_HOME}/bin:${PATH}"
    export JAVA_CMD="${java_bin}"
    return 0
}

bactflow_try_java_home() {
    local home="$1"

    [ -n "${home}" ] || return 1
    bactflow_use_java_bin "${home}/bin/java" && return 0
    bactflow_use_java_bin "${home}/lib/jvm/bin/java" && return 0
    return 1
}

bactflow_pick_nextflow_java() {
    local candidate

    unset JAVA_CMD

    if [ -n "${JAVA_HOME:-}" ]; then
        bactflow_try_java_home "${JAVA_HOME}" && return 0
        bactflow_use_java_bin "${JAVA_HOME}/bin/java" && return 0
    fi

    if [ -n "${BACTFLOW_JAVA_HOME:-}" ]; then
        bactflow_try_java_home "${BACTFLOW_JAVA_HOME}" && return 0
    fi

    for candidate in \
        /usr/lib/jvm/bactflow-java/bin/java \
        /usr/lib/jvm/java-21-openjdk-amd64/bin/java \
        /usr/lib/jvm/java-21-openjdk-arm64/bin/java \
        /usr/lib/jvm/java-17-openjdk-amd64/bin/java \
        /usr/lib/jvm/java-17-openjdk-arm64/bin/java \
        /usr/lib/jvm/default-java/bin/java \
        "${HOME}/.sdkman/candidates/java/current/bin/java"
    do
        bactflow_use_java_bin "${candidate}" && return 0
    done

    for candidate in /usr/lib/jvm/*/bin/java; do
        [ -e "${candidate}" ] || continue
        bactflow_java_in_conda_env "${candidate}" && continue
        bactflow_use_java_bin "${candidate}" && return 0
    done

    if [ -n "${CONDA_PREFIX:-}" ]; then
        bactflow_try_java_home "${CONDA_PREFIX}" && return 0
    fi

    if [ -d /opt/conda/envs/bactflow ]; then
        bactflow_try_java_home /opt/conda/envs/bactflow && return 0
    fi

    return 1
}

bactflow_ensure_nextflow_java() {
    if bactflow_pick_nextflow_java; then
        return 0
    fi

    echo "ERROR: Nextflow needs Java 17-24." >&2
    echo "  - Recreate the bactflow env: mamba env create -f config.yml (pins openjdk=21)" >&2
    echo "  - Or set BACTFLOW_JAVA_HOME to a JDK 17-24 install" >&2
    echo "  - Then run: bash scripts/install_bactflow_hooks.sh" >&2
    return 1
}

bactflow_resolve_nextflow_bin() {
    local candidate

    for candidate in \
        "${BACTFLOW_NEXTFLOW_BIN:-}" \
        "${CONDA_PREFIX:+$CONDA_PREFIX/bin/nextflow}" \
        /opt/conda/envs/bactflow/bin/nextflow
    do
        [ -n "${candidate}" ] && [ -x "${candidate}" ] || continue
        echo "${candidate}"
        return 0
    done

    candidate="$(command -v nextflow 2>/dev/null || true)"
    [ -n "${candidate}" ] && [ -x "${candidate}" ] && echo "${candidate}" && return 0

    return 1
}

bactflow_nextflow_major() {
    local nf_bin="$1"
    NXF_VER="${BACTFLOW_NXF_PIN}" "${nf_bin}" -version 2>&1 \
        | sed -n 's/.*version \([0-9]\+\)\..*/\1/p' \
        | head -n1
}

bactflow_ensure_nextflow_version() {
    local nf_bin major

    nf_bin="$(bactflow_resolve_nextflow_bin)" || {
        echo "ERROR: nextflow not found. Install: mamba install -n bactflow -c bioconda 'nextflow=${BACTFLOW_NXF_PIN}'" >&2
        return 1
    }

    major="$(bactflow_nextflow_major "${nf_bin}")"
    if [ -z "${major}" ] || [ "${major}" -ge "${BACTFLOW_NXF_MAX_MAJOR}" ]; then
        echo "ERROR: BactFlow requires Nextflow < 26 (use ${BACTFLOW_NXF_PIN})." >&2
        echo "  Found: ${nf_bin} (major version: ${major:-unknown})" >&2
        echo "  Fix: mamba install -n bactflow -c bioconda 'nextflow=${BACTFLOW_NXF_PIN}'" >&2
        return 1
    fi

    export BACTFLOW_NEXTFLOW_BIN="${nf_bin}"
    export NXF_VER="${BACTFLOW_NXF_PIN}"
    export PATH="$(dirname "${nf_bin}"):${PATH}"
    return 0
}

bactflow_prepare_nextflow() {
    bactflow_ensure_nextflow_java && bactflow_ensure_nextflow_version
}

if [ -n "${CONDA_PREFIX:-}" ]; then
    _bactflow_hook="${BASH_SOURCE[0]:-$0}"
    case "${_bactflow_hook}" in
        */etc/conda/activate.d/*)
            bactflow_prepare_nextflow || \
                echo "WARNING: bactflow activated but Nextflow/Java setup failed." >&2
            ;;
    esac
fi
