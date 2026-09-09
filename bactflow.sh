#!/bin/bash
# Run BactFlow UI modules in locally built, slim Docker images.

set -euo pipefail

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat << EOF
${GREEN}BactFlow Module Runner${NC}
Pull the latest Docker Hub image for each module (or build locally with --local).

${YELLOW}USAGE:${NC}
    $0 <module> <work_dir> [options]

${YELLOW}MODULES:${NC}
    preassem    Flask UI on port 5000 (seqkit / plots)
    assem       Flask UI on port 5002 (Flye / SPAdes / Unicycler / QUAST)
    postassem   Flask UI on port 5001 (uses the existing post-assembly image)

${YELLOW}REQUIRED:${NC}
    work_dir    Absolute working directory mounted into the container.
                FASTQ folders and --out_dir must live under this path.

${YELLOW}OPTIONS:${NC}
    --no-browser        Don't automatically open a browser
    --install-browser-hook
                        Install/start the host browser hook (for raw docker run)
    --port PORT         Host port (default: 5000 / 5002 / 5001)
    --cpus N            CPU limit (default: auto for preassem, 10 for assem/postassem)
    --memory SIZE       Memory limit (default: auto for preassem, 16g for assem/postassem)
    --pull              Re-pull even if the image is already present (default already refreshes)
    --local             Build/run from the local Dockerfile instead of Docker Hub
    --rebuild           Rebuild the local image (implies --local)
    --tag TAG           Use a specific Hub tag (e.g. v1.0) instead of auto-latest
    --help              Show this help message

${YELLOW}EXAMPLES:${NC}
    $0 preassem /home/user/work_dir
    $0 assem /home/user/work_dir --cpus 16 --memory 32g --port 5002
    $0 postassem /home/user/work_dir --pull
    $0 assem /home/user/work_dir --local --rebuild

EOF
    exit 1
}

ensure_browser_hook() {
    # Tiny host daemon so raw `docker run` containers can open the desktop browser.
    local hook="${SCRIPT_DIR}/scripts/bactflow_browser_hook.sh"
    if [[ ! -f "${hook}" ]]; then
        return 0
    fi
    chmod +x "${hook}" 2>/dev/null || true
    "${hook}" --install >/dev/null 2>&1 || true
}

open_browser() {
    local url=$1
    echo -e "${GREEN}Opening browser to: $url${NC}"

    # Restore GUI session vars when this script is launched from a limited shell.
    if [[ -z "${DISPLAY:-}" ]]; then
        if [[ -S /tmp/.X11-unix/X0 ]]; then
            export DISPLAY=:0
        elif [[ -S /tmp/.X11-unix/X1 ]]; then
            export DISPLAY=:1
        fi
    fi
    if [[ -z "${XDG_RUNTIME_DIR:-}" && -d "/run/user/$(id -u)" ]]; then
        export XDG_RUNTIME_DIR="/run/user/$(id -u)"
    fi
    if [[ -z "${DBUS_SESSION_BUS_ADDRESS:-}" && -n "${XDG_RUNTIME_DIR:-}" && -S "${XDG_RUNTIME_DIR}/bus" ]]; then
        export DBUS_SESSION_BUS_ADDRESS="unix:path=${XDG_RUNTIME_DIR}/bus"
    fi

    # Prefer the same mechanism as local Flask runs (works well on Linux desktops).
    if command -v python3 >/dev/null 2>&1; then
        if python3 - "$url" <<'PY' >/dev/null 2>&1 &
import sys, webbrowser
webbrowser.open(sys.argv[1], new=2, autoraise=True)
PY
        then
            return 0
        fi
    fi

    local opener
    for opener in xdg-open gio open gnome-open sensible-browser firefox google-chrome chromium-browser chromium microsoft-edge; do
        if ! command -v "${opener}" >/dev/null 2>&1; then
            continue
        fi
        if [[ "${opener}" == "gio" ]]; then
            nohup gio open "${url}" >/dev/null 2>&1 &
        else
            nohup "${opener}" "${url}" >/dev/null 2>&1 &
        fi
        return 0
    done

    echo -e "${YELLOW}Could not auto-open a browser. Please visit: ${url}${NC}"
    return 1
}

detect_resources() {
    if [[ "${OSTYPE}" == "linux-gnu"* ]]; then
        AVAILABLE_CPU="$(nproc)"
        AVAILABLE_MEM_G="$(free -g | awk '/^Mem:/{print $2}')"
    elif [[ "${OSTYPE}" == "darwin"* ]]; then
        AVAILABLE_CPU="$(sysctl -n hw.ncpu)"
        AVAILABLE_MEM_G="$(( $(sysctl -n hw.memsize) / 1073741824 ))"
    else
        AVAILABLE_CPU=4
        AVAILABLE_MEM_G=8
    fi
    AVAILABLE_MEM_G="${AVAILABLE_MEM_G:-8}"
}

port_in_use() {
    local port="$1"
    if command -v ss >/dev/null 2>&1; then
        ss -ltn | awk '{print $4}' | grep -Eq "[:.]${port}$"
    elif command -v lsof >/dev/null 2>&1; then
        lsof -iTCP:"${port}" -sTCP:LISTEN >/dev/null 2>&1
    else
        return 1
    fi
}

ensure_docker() {
    if ! command -v docker >/dev/null 2>&1; then
        echo -e "${RED}Error:${NC} Docker is not installed or not on PATH."
        exit 1
    fi
    if ! docker info >/dev/null 2>&1; then
        echo -e "${RED}Error:${NC} Docker daemon is not reachable. Start Docker and add your user to the docker group."
        exit 1
    fi
}

image_exists() {
    docker image inspect "$1" >/dev/null 2>&1
}

build_image() {
    local image="$1"
    local dockerfile="$2"
    local context="$3"
    echo -e "${BLUE}Building ${image} (this is only needed once unless you pass --rebuild)...${NC}"
    docker build \
        -t "${image}" \
        -f "${dockerfile}" \
        "${context}"
}

# Resolve newest version tag on Docker Hub for farhadm1990/<repo>
resolve_latest_hub_tag() {
    local repo="$1"
    local preferred="$2"
    local api="https://hub.docker.com/v2/repositories/farhadm1990/${repo}/tags?page_size=100&ordering=-last_updated"
    local payload tags
    payload="$(curl -fsSL --max-time 20 "${api}" 2>/dev/null || true)"
    tags="$(PAYLOAD="${payload}" python3 - <<'PY'
import json, os, re, sys
raw = os.environ.get("PAYLOAD") or ""
if not raw.strip():
    sys.exit(1)
try:
    data = json.loads(raw)
except Exception:
    sys.exit(1)
tags = []
for row in data.get("results") or []:
    name = row.get("name") or ""
    if name in ("latest", "local"):
        continue
    if re.search(r"\d", name):
        tags.append(name)
if not tags:
    sys.exit(1)

def key(tag):
    m = re.search(r"(\d+(?:\.\d+)*)", tag)
    if not m:
        return (0,)
    return tuple(int(x) for x in m.group(1).split("."))

tags.sort(key=key, reverse=True)
print(tags[0])
PY
)" || true

    if [[ -n "${tags}" ]]; then
        echo "${tags}"
        return 0
    fi
    if [[ -n "${preferred}" ]]; then
        echo "${preferred}"
        return 0
    fi
    echo "latest"
}

pull_hub_image() {
    local image="$1"
    echo -e "${BLUE}Pulling ${image} ...${NC}"
    if docker pull "${image}"; then
        return 0
    fi
    echo -e "${YELLOW}Failed to pull ${image}${NC}"
    return 1
}

resolve_module_image() {
    # Sets IMAGE for the selected module (Hub latest by default).
    local hub_repo="$1"
    local local_image="$2"
    local dockerfile="$3"
    local context="$4"
    local fallback_tag="${5:-}"

    if [[ "${USE_LOCAL}" == true ]]; then
        IMAGE="${local_image}"
        if [[ "${REBUILD}" == true ]] || ! image_exists "${IMAGE}"; then
            if [[ -z "${dockerfile}" || ! -f "${dockerfile}" ]]; then
                echo -e "${RED}Error:${NC} Local Dockerfile missing for ${MODULE}"
                exit 1
            fi
            build_image "${IMAGE}" "${dockerfile}" "${context}"
        else
            echo -e "${BLUE}Using local image ${IMAGE}${NC}"
        fi
        return 0
    fi

    local tag="${IMAGE_TAG}"
    if [[ -z "${tag}" ]]; then
        echo -e "${BLUE}Resolving latest Docker Hub tag for farhadm1990/${hub_repo} ...${NC}"
        tag="$(resolve_latest_hub_tag "${hub_repo}" "${fallback_tag}")"
        echo -e "${GREEN}Latest tag: ${tag}${NC}"
    else
        echo -e "${BLUE}Using requested Hub tag: ${tag}${NC}"
    fi

    IMAGE="farhadm1990/${hub_repo}:${tag}"

    # Always pull the resolved tag so digest stays current when Hub updates the same tag.
    if ! pull_hub_image "${IMAGE}"; then
        if image_exists "${IMAGE}"; then
            echo -e "${YELLOW}Pull failed; using already-cached ${IMAGE}${NC}"
            return 0
        fi
        if [[ -n "${dockerfile}" && -f "${dockerfile}" ]]; then
            echo -e "${YELLOW}Hub pull failed; falling back to local build ${local_image}${NC}"
            IMAGE="${local_image}"
            build_image "${IMAGE}" "${dockerfile}" "${context}"
            return 0
        fi
        echo -e "${RED}Error:${NC} Could not pull ${IMAGE} and no local Dockerfile fallback is available."
        exit 1
    fi
}

OPEN_BROWSER=true
HOST_PORT=""
MODULE=""
WORK_DIR=""
USER_CPUS=""
USER_MEMORY=""
REBUILD=false
USE_LOCAL=false
FORCE_PULL=false
IMAGE_TAG=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --no-browser)
            OPEN_BROWSER=false
            shift
            ;;
        --install-browser-hook)
            ensure_browser_hook
            echo -e "${GREEN}Browser hook installed/started (port ${BACTFLOW_BROWSER_HOOK_PORT:-19264}).${NC}"
            echo "For raw docker run, map -p HOST:CONTAINER and optionally:"
            echo "  -e BACTFLOW_HOST_PORT=HOST -e BACTFLOW_PUBLIC_URL=http://127.0.0.1:HOST/"
            exit 0
            ;;
        --rebuild)
            REBUILD=true
            USE_LOCAL=true
            shift
            ;;
        --local)
            USE_LOCAL=true
            shift
            ;;
        --pull)
            FORCE_PULL=true
            shift
            ;;
        --tag)
            IMAGE_TAG="$2"
            shift 2
            ;;
        --port)
            HOST_PORT="$2"
            shift 2
            ;;
        --cpus)
            USER_CPUS="$2"
            shift 2
            ;;
        --memory)
            USER_MEMORY="$2"
            shift 2
            ;;
        --help|-h)
            usage
            ;;
        -*)
            echo -e "${RED}Error:${NC} Unknown option $1"
            usage
            ;;
        *)
            if [[ -z "${MODULE}" ]]; then
                MODULE="$1"
            elif [[ -z "${WORK_DIR}" ]]; then
                WORK_DIR="$1"
            else
                echo -e "${RED}Error:${NC} Too many arguments"
                usage
            fi
            shift
            ;;
    esac
done

if [[ -z "${MODULE}" || -z "${WORK_DIR}" ]]; then
    echo -e "${RED}Error:${NC} Module and work directory are required"
    usage
fi

ensure_docker
detect_resources

WORK_DIR="$(realpath "${WORK_DIR}")"
if [[ ! -d "${WORK_DIR}" ]]; then
    echo -e "${RED}Error:${NC} Directory '${WORK_DIR}' does not exist"
    exit 1
fi

DOCKERFILE=""
HUB_REPO=""
LOCAL_IMAGE=""
FALLBACK_TAG=""
case "${MODULE}" in
    preassem)
        HUB_REPO="bactflow_preassem"
        LOCAL_IMAGE="bactflow/preassem:local"
        DOCKERFILE="${SCRIPT_DIR}/UI/pre_assem_app/Dockerfile"
        BUILD_CONTEXT="${SCRIPT_DIR}/UI/pre_assem_app"
        FALLBACK_TAG="v1.0"
        DEFAULT_CPUS=""
        DEFAULT_MEMORY=""
        CONTAINER_PORT="5000"
        DEFAULT_HOST_PORT="5000"
        ;;
    assem)
        HUB_REPO="bactflow_assem"
        LOCAL_IMAGE="bactflow/assem:local"
        DOCKERFILE="${SCRIPT_DIR}/UI/assem_app/Dockerfile"
        BUILD_CONTEXT="${SCRIPT_DIR}/UI/assem_app"
        FALLBACK_TAG="v1.0"
        DEFAULT_CPUS="10"
        DEFAULT_MEMORY="16g"
        CONTAINER_PORT="5002"
        DEFAULT_HOST_PORT="5002"
        ;;
    postassem)
        HUB_REPO="bactflow_postassem"
        LOCAL_IMAGE="bactflow/postassem:local"
        DOCKERFILE="${SCRIPT_DIR}/UI/post_assem_app/Dockerfile"
        BUILD_CONTEXT="${SCRIPT_DIR}/UI/post_assem_app"
        FALLBACK_TAG="v0.01"
        DEFAULT_CPUS="10"
        DEFAULT_MEMORY="16g"
        CONTAINER_PORT="5001"
        DEFAULT_HOST_PORT="5001"
        ;;
    *)
        echo -e "${RED}Error:${NC} Unknown module '${MODULE}'"
        usage
        ;;
esac

if [[ -n "${USER_CPUS}" ]]; then
    CPUS="--cpus=${USER_CPUS}"
    CPU_VALUE="${USER_CPUS}"
elif [[ -n "${DEFAULT_CPUS}" ]]; then
    CPUS="--cpus=${DEFAULT_CPUS}"
    CPU_VALUE="${DEFAULT_CPUS}"
else
    CPUS=""
    CPU_VALUE="auto (${AVAILABLE_CPU} cores)"
fi

if [[ -n "${USER_MEMORY}" ]]; then
    MEM_VALUE="${USER_MEMORY}"
elif [[ -n "${DEFAULT_MEMORY}" ]]; then
    MEM_VALUE="${DEFAULT_MEMORY}"
else
    MEM_VALUE=""
fi

if [[ -n "${MEM_VALUE}" ]]; then
    REQ_G="$(echo "${MEM_VALUE}" | sed -E 's/[Gg]$//')"
    if [[ "${REQ_G}" =~ ^[0-9]+$ ]] && [[ "${AVAILABLE_MEM_G}" =~ ^[0-9]+$ ]] && (( REQ_G > AVAILABLE_MEM_G )); then
        CAPPED="$(( AVAILABLE_MEM_G > 2 ? AVAILABLE_MEM_G - 1 : AVAILABLE_MEM_G ))"
        echo -e "${YELLOW}Requested ${MEM_VALUE} RAM but this host has ${AVAILABLE_MEM_G}g. Using ${CAPPED}g.${NC}"
        MEM_VALUE="${CAPPED}g"
    fi
    MEMORY="--memory=${MEM_VALUE}"
else
    MEMORY=""
    MEM_VALUE="auto (${AVAILABLE_MEM_G}g)"
fi

if [[ -z "${HOST_PORT}" ]]; then
    HOST_PORT="${DEFAULT_HOST_PORT}"
fi

if port_in_use "${HOST_PORT}"; then
    echo -e "${RED}Error:${NC} Port ${HOST_PORT} is already in use. Stop the other process or pass --port."
    exit 1
fi

resolve_module_image \
    "${HUB_REPO}" \
    "${LOCAL_IMAGE}" \
    "${DOCKERFILE}" \
    "${BUILD_CONTEXT:-$(dirname "${DOCKERFILE}")}" \
    "${FALLBACK_TAG}"

CONTAINER_NAME="bactflow_${MODULE}_$$"
LOG_FILE="/tmp/${CONTAINER_NAME}.log"
: > "${LOG_FILE}"

echo -e "\n${GREEN}=== Running ${MODULE} ===${NC}"
echo "Work directory: ${WORK_DIR}"
echo "Image: ${IMAGE}"
echo "Resources: CPU=${CPU_VALUE}, Memory=${MEM_VALUE}"
echo "Port mapping: ${HOST_PORT} → ${CONTAINER_PORT}"
echo -e "${CYAN}Access URL: http://127.0.0.1:${HOST_PORT}${NC}"
echo -e "${YELLOW}Press Ctrl+C to stop${NC}\n"

# Keep a host browser hook alive for this session and for raw `docker run` users.
if [[ "${OPEN_BROWSER}" == true ]]; then
    ensure_browser_hook
fi

cleanup() {
    echo -e "\n${YELLOW}Stopping container...${NC}"
    docker stop "${CONTAINER_NAME}" >/dev/null 2>&1 || true
    rm -f "${LOG_FILE}"
}
trap cleanup INT TERM EXIT

# host.docker.internal lets the container reach the browser hook on the host.
DOCKER_EXTRA_HOST=()
if docker run --help 2>/dev/null | grep -q -- '--add-host'; then
    DOCKER_EXTRA_HOST=(--add-host=host.docker.internal:host-gateway)
fi

docker run --rm \
    --name "${CONTAINER_NAME}" \
    --init \
    ${CPUS} \
    ${MEMORY} \
    "${DOCKER_EXTRA_HOST[@]}" \
    -p "${HOST_PORT}:${CONTAINER_PORT}" \
    -v "${WORK_DIR}:${WORK_DIR}" \
    -w "${WORK_DIR}" \
    -e BACTFLOW_IN_DOCKER=1 \
    -e BACTFLOW_NO_BROWSER=1 \
    -e BACTFLOW_MODULE="${MODULE}" \
    -e BACTFLOW_HOST_PORT="${HOST_PORT}" \
    -e BACTFLOW_PUBLIC_URL="http://127.0.0.1:${HOST_PORT}/" \
    -e BACTFLOW_DOCKER_MEMORY="${MEM_VALUE}" \
    -e HOME="${WORK_DIR}" \
    -e NXF_HOME="${WORK_DIR}/.nextflow" \
    "${IMAGE}" > "${LOG_FILE}" 2>&1 &
DOCKER_PID=$!

echo -e "${BLUE}Waiting for the UI to become ready...${NC}"
READY=false
for i in $(seq 1 60); do
    if curl -fsS -o /dev/null "http://127.0.0.1:${HOST_PORT}/"; then
        READY=true
        break
    fi
    if ! kill -0 "${DOCKER_PID}" 2>/dev/null; then
        echo -e "${RED}Container exited before the UI started.${NC}"
        cat "${LOG_FILE}"
        exit 1
    fi
    sleep 1
done

if [[ "${READY}" != true ]]; then
    echo -e "${RED}Timed out waiting for http://127.0.0.1:${HOST_PORT}/${NC}"
    echo "Last container logs:"
    tail -n 80 "${LOG_FILE}" || true
    exit 1
fi

echo -e "${GREEN}UI is ready at http://127.0.0.1:${HOST_PORT}${NC}"
if [[ "${OPEN_BROWSER}" == true ]]; then
    # Slight delay so the first HTTP response is fully settled before the tab opens.
    sleep 0.5
    open_browser "http://127.0.0.1:${HOST_PORT}/" || true
else
    echo -e "${YELLOW}Browser auto-open disabled (--no-browser). Open: http://127.0.0.1:${HOST_PORT}/${NC}"
fi

echo -e "\n${BLUE}Container logs:${NC}\n"
tail -f "${LOG_FILE}" &
TAIL_PID=$!
wait "${DOCKER_PID}"
kill "${TAIL_PID}" >/dev/null 2>&1 || true
echo -e "\n${GREEN}=== ${MODULE} stopped ===${NC}"
