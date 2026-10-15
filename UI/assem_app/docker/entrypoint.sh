#!/usr/bin/env bash
set -euo pipefail

export BACTFLOW_IN_DOCKER=1
export BACTFLOW_MODULE="${BACTFLOW_MODULE:-ui}"
export BACTFLOW_HOST_PORT="${BACTFLOW_HOST_PORT:-5000}"
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

# Prefer an explicit public URL (bactflow.sh / compose set 127.0.0.1 with -p).
# For a raw `docker run` without -p, only the container bridge IP is reachable
# from the host — pick that so the browser opens one working link.
if [ -z "${BACTFLOW_PUBLIC_URL:-}" ]; then
    port="${BACTFLOW_HOST_PORT:-5000}"
    ip=""
    if command -v hostname >/dev/null 2>&1; then
        ip="$(hostname -I 2>/dev/null | awk '{print $1}')"
    fi
    if [ -z "${ip}" ] && [ -r /proc/net/fib_trie ]; then
        ip="$(awk '/32 host/ { print f } {f=$2}' /proc/net/fib_trie 2>/dev/null \
            | grep -vE '^(127\.|0\.0\.0\.0$)' | head -n1 || true)"
    fi
    if [ -n "${ip}" ] && [ "${ip}" != "127.0.0.1" ]; then
        export BACTFLOW_PUBLIC_URL="http://${ip}:${port}/"
    else
        export BACTFLOW_PUBLIC_URL="http://127.0.0.1:${port}/"
    fi
fi
case "${BACTFLOW_PUBLIC_URL}" in
    */) ;;
    *) export BACTFLOW_PUBLIC_URL="${BACTFLOW_PUBLIC_URL}/" ;;
esac

echo ""
echo "=== BactFlow UI ==="
echo "Open this URL: ${BACTFLOW_PUBLIC_URL}"
echo "(Ignore other 'Running on http://…' lines from Flask — use the URL above.)"
echo ""

if [ "$#" -eq 0 ]; then
    echo "No command given to Docker entrypoint" >&2
    exit 1
fi

exec "$@"
