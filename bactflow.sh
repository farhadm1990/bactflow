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
Build (if needed) and run slim pre-assembly / assembly Docker images.

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
    --port PORT         Host port (default: 5000 / 5002 / 5001)
    --cpus N            CPU limit (default: auto for preassem, 10 for assem/postassem)
    --memory SIZE       Memory limit (default: auto for preassem, 16g for assem/postassem)
    --rebuild           Rebuild the image even if it already exists
    --help              Show this help message

${YELLOW}EXAMPLES:${NC}
    $0 preassem /home/user/work_dir
    $0 assem /home/user/work_dir --cpus 16 --memory 32g --port 5002
    $0 preassem /home/user/work_dir --rebuild

EOF
    exit 1
}

open_browser() {
    local url=$1
    echo -e "${GREEN}Opening browser to: $url${NC}"
    if [[ "${OSTYPE}" == "linux-gnu"* ]]; then
        (xdg-open "$url" || gnome-open "$url" || firefox "$url" || true) >/dev/null 2>&1 &
    elif [[ "${OSTYPE}" == "darwin"* ]]; then
        open "$url" >/dev/null 2>&1 &
    else
        echo -e "${YELLOW}Please open: $url${NC}"
    fi
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

OPEN_BROWSER=true
HOST_PORT=""
MODULE=""
WORK_DIR=""
USER_CPUS=""
USER_MEMORY=""
REBUILD=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --no-browser)
            OPEN_BROWSER=false
            shift
            ;;
        --rebuild)
            REBUILD=true
            shift
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
case "${MODULE}" in
    preassem)
        IMAGE="bactflow/preassem:local"
        DOCKERFILE="${SCRIPT_DIR}/UI/pre_assem_app/Dockerfile"
        BUILD_CONTEXT="${SCRIPT_DIR}/UI/pre_assem_app"
        DEFAULT_CPUS=""
        DEFAULT_MEMORY=""
        CONTAINER_PORT="5000"
        DEFAULT_HOST_PORT="5000"
        ;;
    assem)
        IMAGE="bactflow/assem:local"
        DOCKERFILE="${SCRIPT_DIR}/UI/assem_app/Dockerfile"
        BUILD_CONTEXT="${SCRIPT_DIR}/UI/assem_app"
        DEFAULT_CPUS="10"
        DEFAULT_MEMORY="16g"
        CONTAINER_PORT="5002"
        DEFAULT_HOST_PORT="5002"
        ;;
    postassem)
        IMAGE="bactflow/postassem:local"
        DOCKERFILE="${SCRIPT_DIR}/UI/post_assem_app/Dockerfile"
        BUILD_CONTEXT="${SCRIPT_DIR}/UI/post_assem_app"
        FALLBACK_IMAGE="farhadm1990/bactflow_postassem:v0.01"
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

if [[ -n "${DOCKERFILE}" && -f "${DOCKERFILE}" ]]; then
    if [[ "${REBUILD}" == true ]] || ! image_exists "${IMAGE}"; then
        build_image "${IMAGE}" "${DOCKERFILE}" "${BUILD_CONTEXT:-$(dirname "${DOCKERFILE}")}"
    fi
elif [[ -n "${FALLBACK_IMAGE:-}" ]]; then
    IMAGE="${FALLBACK_IMAGE}"
    docker pull "${IMAGE}" >/dev/null || true
else
    echo -e "${RED}Error:${NC} No Dockerfile for ${MODULE}"
    exit 1
fi

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

cleanup() {
    echo -e "\n${YELLOW}Stopping container...${NC}"
    docker stop "${CONTAINER_NAME}" >/dev/null 2>&1 || true
    rm -f "${LOG_FILE}"
}
trap cleanup INT TERM EXIT

docker run --rm \
    --name "${CONTAINER_NAME}" \
    --init \
    ${CPUS} \
    ${MEMORY} \
    -p "${HOST_PORT}:${CONTAINER_PORT}" \
    -v "${WORK_DIR}:${WORK_DIR}" \
    -w "${WORK_DIR}" \
    -e BACTFLOW_IN_DOCKER=1 \
    -e BACTFLOW_NO_BROWSER=1 \
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
    open_browser "http://127.0.0.1:${HOST_PORT}"
fi

echo -e "\n${BLUE}Container logs:${NC}\n"
tail -f "${LOG_FILE}" &
TAIL_PID=$!
wait "${DOCKER_PID}"
kill "${TAIL_PID}" >/dev/null 2>&1 || true
echo -e "\n${GREEN}=== ${MODULE} stopped ===${NC}"
