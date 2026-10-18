#!/bin/bash
# Launch the unified BactFlow UI with Docker, Apptainer/Singularity, or Conda.

set -euo pipefail

# $'…' stores a real ESC. Single-quoted '\033' printed the codes as text in --help.
if [[ -t 1 ]]; then
    GREEN=$'\033[32m'
    RED=$'\033[31m'
    YELLOW=$'\033[33m'
    BLUE=$'\033[34m'
    CYAN=$'\033[36m'
    NC=$'\033[0m'
else
    GREEN="" RED="" YELLOW="" BLUE="" CYAN="" NC=""
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat << EOF
${GREEN}BactFlow${NC}
Start the unified UI (Home, QC, Assembly, Post-assembly, Wiki).
Do not call docker / singularity / python yourself.

${YELLOW}USAGE:${NC}
    $0 run <work_dir> [options]

${YELLOW}REQUIRED:${NC}
    work_dir    Working directory (mounted into containers).
                FASTQ folders and --out_dir should live under this path.

${YELLOW}OPTIONS:${NC}
    --service NAME      Runtime: docker (default), apptainer, or conda
    --no-browser        Don't automatically open a browser
    --port PORT         Host port (default: 5000)
    --cpus N            CPU limit (Docker / Apptainer when supported)
    --memory SIZE       Memory limit (Docker; Apptainer when supported)
    --pull              Re-pull the Hub image / rebuild the converted Apptainer image
    --local             Build/run from the local Dockerfile (docker / apptainer)
    --rebuild           Rebuild the local Docker image, Apptainer image, or Conda env
    --tag TAG           Pin a Hub tag (e.g. v1.0) instead of auto-latest
    --help              Show this help message

    Env: BACTFLOW_PROFILE=full|slim  (local docker build; slim omits geNomad ~half size)

${YELLOW}EXAMPLES:${NC}
    $0 run /home/user/work_dir
    $0 run /home/user/work_dir --service conda
    $0 run /home/user/work_dir --service docker --local --rebuild
    BACTFLOW_PROFILE=slim $0 run /home/user/work_dir --service docker --local --rebuild
    $0 run /home/user/work_dir --service docker --cpus 16 --memory 32g
    $0 run /home/user/work_dir --service apptainer
    docker compose up --build

EOF
    exit 1
}

ensure_browser_hook() {
    # Tiny host daemon so the Docker service can open the desktop browser.
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

normalize_service() {
    local raw="${1,,}"
    case "${raw}" in
        docker|podman) echo "docker" ;;
        apptainer|singularity|aptainer) echo "apptainer" ;;
        conda|mamba) echo "conda" ;;
        *) echo "" ;;
    esac
}

legacy_module_name() {
    case "${1}" in
        pre|preassem|pre-assembly|pre_assem|preassembly|ass|assem|assembly|post|postassem|post-assembly|post_assem|postassembly)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

ensure_docker() {
    if ! command -v docker >/dev/null 2>&1; then
        echo -e "${RED}Error:${NC} Docker is not installed or not on PATH."
        echo "Use --service apptainer or --service conda if you cannot use Docker."
        exit 1
    fi
    if ! docker info >/dev/null 2>&1; then
        echo -e "${RED}Error:${NC} Docker daemon is not reachable. Start Docker and add your user to the docker group."
        echo "Sites without Docker privileges can use --service apptainer or --service conda."
        exit 1
    fi
}

ensure_apptainer() {
    if command -v apptainer >/dev/null 2>&1; then
        APPTAINER_BIN="apptainer"
        return 0
    fi
    if command -v singularity >/dev/null 2>&1; then
        APPTAINER_BIN="singularity"
        return 0
    fi
    echo -e "${RED}Error:${NC} Apptainer/Singularity is not installed or not on PATH."
    echo "Install Apptainer (no docker group needed), or use --service conda."
    exit 1
}

# Singularity 3.0 has no --env; APPTAINERENV_/SINGULARITYENV_ work on old and new.
apptainer_export() {
    local key="$1"
    local val="$2"
    export "${key}=${val}"
    export "SINGULARITYENV_${key}=${val}"
    export "APPTAINERENV_${key}=${val}"
}

apptainer_help_has() {
    local flag="$1"
    grep -qE "(^|[[:space:]])${flag}([[:space:],]|$)" <<<"${APPTAINER_RUN_HELP}"
}

# True if Singularity/Apptainer can actually start this image (file or sandbox dir).
apptainer_image_ok() {
    local img="$1"
    [[ -e "${img}" ]] || return 1
    "${APPTAINER_BIN}" exec "${img}" /bin/true >/dev/null 2>&1
}

apptainer_squashfs_offset() {
    local sif="$1"
    local offset=""
    if command -v python3 >/dev/null 2>&1; then
        offset="$(python3 -c 'import sys; d=open(sys.argv[1],"rb").read(); i=d.find(b"hsqs"); sys.exit(1) if i<0 else print(i)' "${sif}" 2>/dev/null || true)"
    fi
    if [[ -z "${offset}" ]]; then
        offset="$(LC_ALL=C grep -aob $'hsqs' "${sif}" 2>/dev/null | head -1 | cut -d: -f1 || true)"
    fi
    [[ "${offset}" =~ ^[0-9]+$ ]] || return 1
    printf '%s\n' "${offset}"
}

# Old setuid Singularity (3.0) + a new kernel often cannot mount SIF squashfs
# (EINVAL). A sandbox directory uses the same rootfs without that mount.
apptainer_sif_to_sandbox() {
    local sif="$1"
    local dest="$2"
    local offset
    command -v unsquashfs >/dev/null 2>&1 || return 1
    offset="$(apptainer_squashfs_offset "${sif}")" || return 1
    [[ "${offset}" =~ ^[0-9]+$ ]] || return 1
    echo -e "${BLUE}Kernel cannot mount this .sif; extracting a sandbox to ${dest} ...${NC}"
    rm -rf "${dest}"
    unsquashfs -o "${offset}" -f -d "${dest}" "${sif}" || return 1
    [[ -d "${dest}/.singularity.d" ]]
}

find_conda_base() {
    local candidate
    if command -v conda >/dev/null 2>&1; then
        conda info --base 2>/dev/null && return 0
    fi
    for candidate in \
        "${HOME}/miniforge3" \
        "${HOME}/mambaforge" \
        "${HOME}/miniconda3" \
        "${HOME}/anaconda3" \
        "/opt/conda" \
        "/opt/miniconda3"
    do
        if [[ -x "${candidate}/bin/conda" ]]; then
            echo "${candidate}"
            return 0
        fi
    done
    return 1
}

ensure_conda() {
    CONDA_BASE="$(find_conda_base || true)"
    if [[ -z "${CONDA_BASE}" ]]; then
        echo -e "${RED}Error:${NC} conda/mamba was not found."
        echo "Install Miniforge or Miniconda in your home directory (no sudo):"
        echo "  curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
        echo "  bash Miniforge3-Linux-x86_64.sh -b -p \"\$HOME/miniforge3\""
        exit 1
    fi
    # shellcheck disable=SC1091
    set +u
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    # Keep CONDA_EXE as conda. Mamba 2 rejects `shell.posix`, which conda.sh
    # uses for `conda activate`. Use mamba only for env create/update.
    export CONDA_EXE="${CONDA_BASE}/bin/conda"
    export PATH="${CONDA_BASE}/bin:${PATH}"
    set -u
}

conda_env_exists() {
    conda env list 2>/dev/null | awk '{print $1}' | grep -qx "bactflow"
}

image_exists() {
    docker image inspect "$1" >/dev/null 2>&1
}

build_image() {
    local image="$1"
    local dockerfile="$2"
    local context="$3"
    local profile="${BACTFLOW_PROFILE:-full}"
    echo -e "${BLUE}Building ${image} (profile=${profile}; once unless you pass --rebuild)...${NC}"
    docker build \
        --build-arg "BACTFLOW_PROFILE=${profile}" \
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
    IMAGE_TAG="${tag}"
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

# Convert an OCI/Docker image without using Singularity's old Docker API
# client (1.22). Prefer `docker save` → docker-archive:// when Docker works.
# Third argument "sandbox" builds a chroot directory (no squashfs mount).
apptainer_build_from() {
    local dest="$1"
    local src="$2"
    local sandbox="${3:-}"
    local authdir
    authdir="$(mktemp -d "${TMPDIR:-/tmp}/bactflow-apptainer-auth.XXXXXX")"
    printf '%s\n' '{"auths":{}}' > "${authdir}/config.json"
    printf '%s\n' '{"auths":{}}' > "${authdir}/auth.json"
    unset APPTAINER_DOCKER_USERNAME APPTAINER_DOCKER_PASSWORD \
        SINGULARITY_DOCKER_USERNAME SINGULARITY_DOCKER_PASSWORD \
        DOCKER_USERNAME DOCKER_PASSWORD || true
    # Stop Apptainer/Singularity from sending a stale `docker login`.
    export DOCKER_CONFIG="${authdir}"
    export APPTAINER_AUTHFILE="${authdir}/auth.json"
    export SINGULARITY_AUTHFILE="${authdir}/auth.json"

    local args=(build --force)
    if [[ "${sandbox}" == sandbox ]]; then
        args+=(--sandbox)
        rm -rf "${dest}"
    fi
    if "${APPTAINER_BIN}" build --help 2>/dev/null | grep -q -- '--authfile'; then
        args+=(--authfile "${authdir}/auth.json")
    fi
    echo -e "${BLUE}Building ${dest}${NC}"
    echo -e "${BLUE}Source: ${src}${NC}"
    local ok=0
    if "${APPTAINER_BIN}" "${args[@]}" "${dest}" "${src}"; then
        ok=1
    fi
    rm -rf "${authdir}"
    unset APPTAINER_AUTHFILE SINGULARITY_AUTHFILE DOCKER_CONFIG || true
    [[ "${ok}" == 1 ]]
}

# docker pull (modern CLI) then tar → image. Avoids Hub auth inside Singularity
# and avoids docker-daemon:// which needs Docker API ≥ 1.24.
apptainer_build_from_docker_tar() {
    local dest="$1"
    local image="$2"
    local tar_path="$3"
    local sandbox="${4:-}"
    if ! docker pull "${image}"; then
        if ! docker image inspect "${image}" >/dev/null 2>&1; then
            return 1
        fi
        echo -e "${YELLOW}docker pull failed; using the image already on this machine.${NC}"
    fi
    echo -e "${BLUE}docker save → ${tar_path}${NC}"
    docker save -o "${tar_path}" "${image}" || return 1
    local src="docker-archive://${tar_path}"
    # Absolute paths need three slashes: docker-archive:///abs/path.tar
    if [[ "${tar_path}" == /* ]]; then
        src="docker-archive://${tar_path}"
    fi
    apptainer_build_from "${dest}" "${src}" "${sandbox}"
    local ok=$?
    rm -f "${tar_path}"
    return "${ok}"
}

resolve_apptainer_sif() {
    local hub_repo="$1"
    local local_image="$2"
    local dockerfile="$3"
    local context="$4"
    local fallback_tag="${5:-}"
    local sif_dir="${BACTFLOW_SIF_DIR:-${HOME}/.bactflow/sif}"
    mkdir -p "${sif_dir}"

    local tag="${IMAGE_TAG}"
    local src
    local sif_file
    local sandbox
    local mark="${sif_dir}/.prefer_sandbox"
    local want_sandbox=false
    [[ -f "${mark}" ]] && want_sandbox=true

    if [[ "${USE_LOCAL}" == true ]]; then
        ensure_docker
        IMAGE="${local_image}"
        if [[ "${REBUILD}" == true ]] || ! image_exists "${IMAGE}"; then
            if [[ -z "${dockerfile}" || ! -f "${dockerfile}" ]]; then
                echo -e "${RED}Error:${NC} Local Dockerfile missing for ${MODULE}"
                exit 1
            fi
            build_image "${IMAGE}" "${dockerfile}" "${context}"
        fi
        tag="local"
        src="docker-archive"
        sif_file="${sif_dir}/${hub_repo}_local.sif"
        sandbox="${sif_dir}/${hub_repo}_local"
    else
        if [[ -z "${tag}" ]]; then
            echo -e "${BLUE}Resolving latest Docker Hub tag for farhadm1990/${hub_repo} ...${NC}"
            tag="$(resolve_latest_hub_tag "${hub_repo}" "${fallback_tag}")"
            echo -e "${GREEN}Latest tag: ${tag}${NC}"
        else
            echo -e "${BLUE}Using requested Hub tag: ${tag}${NC}"
        fi
        IMAGE_TAG="${tag}"
        IMAGE="farhadm1990/${hub_repo}:${tag}"
        src="docker://${IMAGE}"
        sif_file="${sif_dir}/${hub_repo}_${tag}.sif"
        sandbox="${sif_dir}/${hub_repo}_${tag}"
    fi

    if [[ "${FORCE_PULL}" != true && "${REBUILD}" != true ]]; then
        if [[ -d "${sandbox}" ]] && apptainer_image_ok "${sandbox}"; then
            echo -e "${BLUE}Using existing sandbox ${sandbox}${NC}"
            echo -e "${YELLOW}Pass --pull to rebuild it from ${src}${NC}"
            SIF="${sandbox}"
            return 0
        fi
        if [[ -f "${sif_file}" ]]; then
            if apptainer_image_ok "${sif_file}"; then
                echo -e "${BLUE}Using existing ${sif_file}${NC}"
                echo -e "${YELLOW}Pass --pull to rebuild it from ${src}${NC}"
                SIF="${sif_file}"
                return 0
            fi
            echo -e "${YELLOW}Existing .sif cannot be mounted (typical of Singularity 3.0 on a new kernel).${NC}"
            if apptainer_sif_to_sandbox "${sif_file}" "${sandbox}" && apptainer_image_ok "${sandbox}"; then
                date > "${mark}"
                SIF="${sandbox}"
                return 0
            fi
            want_sandbox=true
        fi
    fi

    local dest="${sif_file}"
    local mode=""
    if [[ "${want_sandbox}" == true ]]; then
        dest="${sandbox}"
        mode="sandbox"
    fi

    local tar_path="${sif_dir}/.${hub_repo}_${tag}.tar"
    if command -v docker >/dev/null 2>&1 && docker info >/dev/null 2>&1; then
        echo -e "${BLUE}Converting Docker image to ${APPTAINER_BIN} via docker save (no Docker daemon API)...${NC}"
        if apptainer_build_from_docker_tar "${dest}" "${IMAGE}" "${tar_path}" "${mode}"; then
            if apptainer_image_ok "${dest}"; then
                SIF="${dest}"
                return 0
            fi
            if [[ -f "${dest}" ]] && apptainer_sif_to_sandbox "${dest}" "${sandbox}" && apptainer_image_ok "${sandbox}"; then
                date > "${mark}"
                SIF="${sandbox}"
                return 0
            fi
            echo -e "${YELLOW}Built image is not runnable; trying sandbox from Docker tar...${NC}"
            if [[ "${mode}" != sandbox ]] && apptainer_build_from_docker_tar "${sandbox}" "${IMAGE}" "${tar_path}" sandbox && apptainer_image_ok "${sandbox}"; then
                date > "${mark}"
                SIF="${sandbox}"
                return 0
            fi
        fi
        echo -e "${YELLOW}docker save conversion failed; trying ${src} ...${NC}"
    fi

    if [[ "${src}" == docker://* ]] && apptainer_build_from "${dest}" "${src}" "${mode}"; then
        if apptainer_image_ok "${dest}"; then
            SIF="${dest}"
            return 0
        fi
        if [[ -f "${dest}" ]] && apptainer_sif_to_sandbox "${dest}" "${sandbox}" && apptainer_image_ok "${sandbox}"; then
            date > "${mark}"
            SIF="${sandbox}"
            return 0
        fi
    fi

    echo -e "${RED}Error:${NC} Could not build a runnable image for ${MODULE}"
    echo "Apptainer/Singularity converts the published Docker Hub image."
    echo "On old Singularity (3.0) a .sif often cannot be mounted; BactFlow then"
    echo "extracts a sandbox directory under ${sif_dir}."
    echo "Need either a working docker CLI, or:"
    echo "  ./bactflow.sh run \"${WORK_DIR}\" --service conda"
    exit 1
}

wait_for_ui() {
    local pid="$1"
    echo -e "${BLUE}Waiting for the UI to become ready...${NC}"
    READY=false
    local i
    for i in $(seq 1 90); do
        if curl -fsS -o /dev/null "http://127.0.0.1:${HOST_PORT}/"; then
            READY=true
            break
        fi
        if ! kill -0 "${pid}" 2>/dev/null; then
            echo -e "${RED}Process exited before the UI started.${NC}"
            cat "${LOG_FILE}"
            exit 1
        fi
        sleep 1
    done

    if [[ "${READY}" != true ]]; then
        echo -e "${RED}Timed out waiting for http://127.0.0.1:${HOST_PORT}/${NC}"
        echo "Last logs:"
        tail -n 80 "${LOG_FILE}" || true
        exit 1
    fi

    echo -e "${GREEN}UI is ready at http://127.0.0.1:${HOST_PORT}${NC}"
    if [[ "${OPEN_BROWSER}" == true ]]; then
        sleep 0.5
        open_browser "http://127.0.0.1:${HOST_PORT}/" || true
    else
        echo -e "${YELLOW}Browser auto-open disabled (--no-browser). Open: http://127.0.0.1:${HOST_PORT}/${NC}"
    fi
}

follow_until_exit() {
    local pid="$1"
    echo -e "\n${BLUE}Logs:${NC}\n"
    tail -f "${LOG_FILE}" &
    TAIL_PID=$!
    wait "${pid}"
    kill "${TAIL_PID}" >/dev/null 2>&1 || true
    echo -e "\n${GREEN}=== BactFlow stopped ===${NC}"
}

# Normalize --memory values: 20 → 20g, 20G/20gb → 20g (Docker treats bare
# integers as bytes, which is never what users mean).
normalize_memory_value() {
    local raw="${1:-}"
    raw="$(echo "${raw}" | tr '[:upper:]' '[:lower:]' | tr -d ' ')"
    if [[ -z "${raw}" || "${raw}" == auto* ]]; then
        echo "${raw}"
        return 0
    fi
    if [[ "${raw}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
        echo "${raw%%.*}g"
        return 0
    fi
    raw="${raw/%gi/g}"
    raw="${raw/%gb/g}"
    raw="${raw/%mi/m}"
    raw="${raw/%mb/m}"
    echo "${raw}"
}

memory_gb_number() {
    local raw
    raw="$(normalize_memory_value "${1:-}")"
    raw="${raw/%g/}"
    raw="${raw/%m/}"
    if [[ "${raw}" =~ ^[0-9]+$ ]]; then
        echo "${raw}"
    else
        echo ""
    fi
}

export_runtime_env() {
    export BACTFLOW_MODULE="${MODULE}"
    export BACTFLOW_HOST_PORT="${HOST_PORT}"
    export BACTFLOW_PUBLIC_URL="http://127.0.0.1:${HOST_PORT}/"
    export BACTFLOW_NO_BROWSER=1
    export BACTFLOW_SERVICE="${SERVICE}"
    export HOME="${WORK_DIR}"
    export NXF_HOME="${WORK_DIR}/.nextflow"
    if [[ -n "${MEM_VALUE}" && "${MEM_VALUE}" != auto* ]]; then
        export BACTFLOW_MEMORY="${MEM_VALUE}"
        export BACTFLOW_DOCKER_MEMORY="${MEM_VALUE}"
        local mem_gb
        mem_gb="$(memory_gb_number "${MEM_VALUE}")"
        if [[ -n "${mem_gb}" ]]; then
            export BACTFLOW_MEMORY_GB="${mem_gb}"
        fi
    fi
    if [[ "${CPU_VALUE}" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
        export BACTFLOW_CPUS="${CPU_VALUE}"
        export NXF_OPTS="${NXF_OPTS:-}"
    fi
}

run_docker() {
    ensure_docker
    resolve_module_image \
        "${HUB_REPO}" \
        "${LOCAL_IMAGE}" \
        "${DOCKERFILE}" \
        "${BUILD_CONTEXT}" \
        "${FALLBACK_TAG}"

    CONTAINER_NAME="bactflow_${MODULE}_$$"
    LOG_FILE="/tmp/${CONTAINER_NAME}.log"
    : > "${LOG_FILE}"

    echo -e "\n${GREEN}=== Running BactFlow (docker) ===${NC}"
    echo "Work directory: ${WORK_DIR}"
    echo "Image: ${IMAGE}"
    echo "Resources: CPU=${CPU_VALUE}, Memory=${MEM_VALUE}"
    echo "Port mapping: ${HOST_PORT} → ${CONTAINER_PORT}"
    echo -e "${CYAN}Access URL: http://127.0.0.1:${HOST_PORT}${NC}"
    echo -e "${YELLOW}Press Ctrl+C to stop${NC}\n"

    if [[ "${OPEN_BROWSER}" == true ]]; then
        ensure_browser_hook
    fi

    cleanup() {
        echo -e "\n${YELLOW}Stopping container...${NC}"
        docker stop "${CONTAINER_NAME}" >/dev/null 2>&1 || true
        rm -f "${LOG_FILE}"
    }
    trap cleanup INT TERM EXIT

    local extra_host=()
    if docker run --help 2>/dev/null | grep -q -- '--add-host'; then
        extra_host=(--add-host=host.docker.internal:host-gateway)
    fi

    # --memory-swap equal to --memory disables extra swap so RAM cannot
    # quietly grow past the requested cap (Docker default is often 2x).
    local MEMORY_SWAP=()
    if [[ -n "${MEMORY}" ]]; then
        MEMORY_SWAP=(--memory-swap="${MEM_VALUE}")
    fi

    docker run --rm \
        --name "${CONTAINER_NAME}" \
        --init \
        ${CPUS} \
        ${MEMORY} \
        "${MEMORY_SWAP[@]}" \
        "${extra_host[@]}" \
        -p "${HOST_PORT}:${CONTAINER_PORT}" \
        -v "${WORK_DIR}:${WORK_DIR}" \
        -w "${WORK_DIR}" \
        -e BACTFLOW_IN_DOCKER=1 \
        -e BACTFLOW_NO_BROWSER=1 \
        -e BACTFLOW_MODULE="${MODULE}" \
        -e BACTFLOW_HOST_PORT="${CONTAINER_PORT}" \
        -e BACTFLOW_PUBLIC_URL="http://127.0.0.1:${HOST_PORT}/" \
        -e BACTFLOW_DOCKER_MEMORY="${MEM_VALUE}" \
        -e BACTFLOW_MEMORY="${MEM_VALUE}" \
        -e BACTFLOW_MEMORY_GB="$(memory_gb_number "${MEM_VALUE}")" \
        -e BACTFLOW_SERVICE=docker \
        -e HOME="${WORK_DIR}" \
        -e NXF_HOME="${WORK_DIR}/.nextflow" \
        "${IMAGE}" > "${LOG_FILE}" 2>&1 &
    local pid=$!
    wait_for_ui "${pid}"
    follow_until_exit "${pid}"
}

run_apptainer() {
    ensure_apptainer
    resolve_apptainer_sif \
        "${HUB_REPO}" \
        "${LOCAL_IMAGE}" \
        "${DOCKERFILE}" \
        "${BUILD_CONTEXT}" \
        "${FALLBACK_TAG}"

    LOG_FILE="/tmp/bactflow_${MODULE}_apptainer_$$.log"
    : > "${LOG_FILE}"

    echo -e "\n${GREEN}=== Running BactFlow (${APPTAINER_BIN}) ===${NC}"
    echo "Work directory: ${WORK_DIR}"
    echo "Image: ${SIF}"
    echo "Resources: CPU=${CPU_VALUE}, Memory=${MEM_VALUE}"
    echo -e "${CYAN}Access URL: http://127.0.0.1:${HOST_PORT}${NC}"
    echo -e "${YELLOW}Press Ctrl+C to stop${NC}\n"

    local APPTAINER_RUN_HELP extra app_in_sif
    APPTAINER_RUN_HELP="$("${APPTAINER_BIN}" exec --help 2>&1 || true)"
    extra=(--bind "${WORK_DIR}:${WORK_DIR}")
    if apptainer_help_has --pwd; then
        extra+=(--pwd "${WORK_DIR}")
    fi
    if apptainer_help_has --home; then
        extra+=(--home "${WORK_DIR}")
    fi
    if [[ "${CPU_VALUE}" =~ ^[0-9]+([.][0-9]+)?$ ]] && apptainer_help_has --cpus; then
        extra+=(--cpus "${CPU_VALUE}")
    fi
    if [[ -n "${MEM_VALUE}" && "${MEM_VALUE}" != auto* ]] && apptainer_help_has --memory; then
        extra+=(--memory "${MEM_VALUE}")
    fi

    export_runtime_env
    # Host network: the process binds HOST_PORT on the machine (no docker -p).
    apptainer_export BACTFLOW_HOST_PORT "${HOST_PORT}"
    apptainer_export BACTFLOW_IN_DOCKER 1
    apptainer_export BACTFLOW_NO_BROWSER 1
    apptainer_export BACTFLOW_MODULE "${MODULE}"
    apptainer_export BACTFLOW_PUBLIC_URL "http://127.0.0.1:${HOST_PORT}/"
    apptainer_export BACTFLOW_SERVICE apptainer
    apptainer_export HOME "${WORK_DIR}"
    apptainer_export NXF_HOME "${WORK_DIR}/.nextflow"
    if [[ -n "${MEM_VALUE}" && "${MEM_VALUE}" != auto* ]]; then
        apptainer_export BACTFLOW_MEMORY "${MEM_VALUE}"
        apptainer_export BACTFLOW_DOCKER_MEMORY "${MEM_VALUE}"
        apptainer_export BACTFLOW_MEMORY_GB "$(memory_gb_number "${MEM_VALUE}")"
    fi
    # Only new Apptainer/Singularity accept --env KEY=VAL (3.0 prints Usage and exits).
    if apptainer_help_has --env; then
        extra+=(--env "BACTFLOW_IN_DOCKER=1")
        extra+=(--env "BACTFLOW_NO_BROWSER=1")
        extra+=(--env "BACTFLOW_MODULE=${MODULE}")
        extra+=(--env "BACTFLOW_HOST_PORT=${HOST_PORT}")
        extra+=(--env "BACTFLOW_PUBLIC_URL=http://127.0.0.1:${HOST_PORT}/")
        extra+=(--env "BACTFLOW_SERVICE=apptainer")
        extra+=(--env "HOME=${WORK_DIR}")
        extra+=(--env "NXF_HOME=${WORK_DIR}/.nextflow")
        if [[ -n "${MEM_VALUE}" && "${MEM_VALUE}" != auto* ]]; then
            extra+=(--env "BACTFLOW_MEMORY=${MEM_VALUE}")
            extra+=(--env "BACTFLOW_DOCKER_MEMORY=${MEM_VALUE}")
            extra+=(--env "BACTFLOW_MEMORY_GB=$(memory_gb_number "${MEM_VALUE}")")
        fi
    fi

    cleanup() {
        echo -e "\n${YELLOW}Stopping ${APPTAINER_BIN} process...${NC}"
        if [[ -n "${APP_PID:-}" ]]; then
            kill "${APP_PID}" >/dev/null 2>&1 || true
            wait "${APP_PID}" >/dev/null 2>&1 || true
        fi
        rm -f "${LOG_FILE}"
    }
    trap cleanup INT TERM EXIT

    # Image layout (unified UI): /app/UI/bactflow_app/app.py (+ symlink /app/app.py).
    app_in_sif="${APP_IN_IMAGE}"
    # exec + entrypoint: old Singularity Docker conversions often have an empty runscript.
    "${APPTAINER_BIN}" exec \
        "${extra[@]}" \
        "${SIF}" \
        /usr/local/bin/bactflow-entrypoint python "${app_in_sif}" > "${LOG_FILE}" 2>&1 &
    APP_PID=$!
    wait_for_ui "${APP_PID}"
    follow_until_exit "${APP_PID}"
}

run_conda() {
    ensure_conda
    if [[ "${USE_LOCAL}" == true || -n "${IMAGE_TAG}" || "${FORCE_PULL}" == true ]]; then
        echo -e "${YELLOW}Note: --local / --tag / --pull apply to docker and apptainer, not conda.${NC}"
    fi

    if [[ "${REBUILD}" == true ]] || ! conda_env_exists; then
        echo -e "${BLUE}Creating/updating the bactflow Conda environment (first time can take a while)...${NC}"
        bash "${SCRIPT_DIR}/scripts/setup_bactflow_env.sh" "${SCRIPT_DIR}/config.yml" bactflow
        # setup_bactflow_env.sh leaves the shell unactivated
        # shellcheck disable=SC1091
        set +u
        source "${CONDA_BASE}/etc/profile.d/conda.sh"
        set -u
    else
        echo -e "${BLUE}Using existing Conda env 'bactflow'${NC}"
        echo -e "${YELLOW}Pass --rebuild to update it from config.yml${NC}"
    fi

    # conda activate scripts are not `set -u` safe (MKL INTERFACE_LAYER backup).
    set +u
    conda activate bactflow
    set -u

    LOG_FILE="/tmp/bactflow_${MODULE}_conda_$$.log"
    : > "${LOG_FILE}"

    echo -e "\n${GREEN}=== Running BactFlow (conda) ===${NC}"
    echo "Work directory: ${WORK_DIR}"
    echo "Env: ${CONDA_PREFIX:-bactflow}"
    echo "App: ${APP_PY}"
    echo "Resources: CPU=${CPU_VALUE}, Memory=${MEM_VALUE}"
    echo -e "${CYAN}Access URL: http://127.0.0.1:${HOST_PORT}${NC}"
    echo -e "${YELLOW}Press Ctrl+C to stop${NC}\n"

    export_runtime_env
    unset BACTFLOW_IN_DOCKER || true

    cleanup() {
        echo -e "\n${YELLOW}Stopping Conda UI...${NC}"
        if [[ -n "${APP_PID:-}" ]]; then
            kill "${APP_PID}" >/dev/null 2>&1 || true
            wait "${APP_PID}" >/dev/null 2>&1 || true
        fi
        rm -f "${LOG_FILE}"
    }
    trap cleanup INT TERM EXIT

    # Prefer a hard MemoryMax scope for the whole UI+pipeline tree when systemd
    # user sessions are available. Pipeline code also wraps Nextflow the same way.
    local run_cmd=(python "${APP_PY}")
    local mem_gb
    mem_gb="$(memory_gb_number "${MEM_VALUE}")"
    if [[ -n "${mem_gb}" ]] && command -v systemd-run >/dev/null 2>&1; then
        if systemd-run --user --quiet --collect --scope -p MemoryMax=64M -p MemorySwapMax=0 -- true >/dev/null 2>&1; then
            echo -e "${BLUE}Hard RAM cap: ${mem_gb}g (systemd MemoryMax)${NC}"
            run_cmd=(
                systemd-run --user --quiet --collect --scope
                -p "MemoryMax=${mem_gb}G"
                -p MemorySwapMax=0
                -- python "${APP_PY}"
            )
        else
            echo -e "${YELLOW}WARNING: systemd user MemoryMax unavailable; pipeline will try a per-run scope. Prefer --service docker --memory ${mem_gb}g for a hard OS cap.${NC}"
        fi
    fi

    (
        cd "${WORK_DIR}"
        "${run_cmd[@]}"
    ) > "${LOG_FILE}" 2>&1 &
    APP_PID=$!
    wait_for_ui "${APP_PID}"
    follow_until_exit "${APP_PID}"
}

OPEN_BROWSER=true
HOST_PORT=""
MODULE="ui"
WORK_DIR=""
CMD=""
USER_CPUS=""
USER_MEMORY=""
REBUILD=false
USE_LOCAL=false
FORCE_PULL=false
IMAGE_TAG=""
SERVICE="docker"
APPTAINER_BIN=""
SIF=""
APP_PY=""
CONDA_BASE=""
CONDA_EXE="conda"

while [[ $# -gt 0 ]]; do
    case $1 in
        --service)
            SERVICE="$(normalize_service "${2:-}")"
            if [[ -z "${SERVICE}" ]]; then
                echo -e "${RED}Error:${NC} --service must be docker, apptainer, or conda"
                usage
            fi
            shift 2
            ;;
        --service=*)
            SERVICE="$(normalize_service "${1#--service=}")"
            if [[ -z "${SERVICE}" ]]; then
                echo -e "${RED}Error:${NC} --service must be docker, apptainer, or conda"
                usage
            fi
            shift
            ;;
        --no-browser)
            OPEN_BROWSER=false
            shift
            ;;
        --install-browser-hook)
            ensure_browser_hook
            echo -e "${GREEN}Browser hook installed/started (port ${BACTFLOW_BROWSER_HOOK_PORT:-19264}).${NC}"
            echo "You do not need this flag for normal use; ./bactflow.sh starts the hook itself."
            exit 0
            ;;
        --rebuild)
            REBUILD=true
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
            if [[ -z "${CMD}" ]]; then
                if [[ "$1" == "run" || "$1" == "start" || "$1" == "ui" ]]; then
                    CMD="run"
                elif legacy_module_name "$1"; then
                    echo -e "${RED}Error:${NC} Separate module launches were removed."
                    echo "Use: $0 run <work_dir> [--service docker|apptainer|conda]"
                    exit 1
                elif [[ -d "$1" ]]; then
                    CMD="run"
                    WORK_DIR="$1"
                else
                    echo -e "${RED}Error:${NC} Unknown argument '$1'"
                    usage
                fi
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

if [[ "${CMD}" != "run" || -z "${WORK_DIR}" ]]; then
    echo -e "${RED}Error:${NC} Usage: $0 run <work_dir> [--service docker|apptainer|conda]"
    usage
fi
# Historic docker behaviour: --rebuild means rebuild the local image.
if [[ "${SERVICE}" == "docker" && "${REBUILD}" == true ]]; then
    USE_LOCAL=true
fi
if [[ "${SERVICE}" == "conda" ]]; then
    USE_LOCAL=false
fi

detect_resources

if [[ ! -d "${WORK_DIR}" ]]; then
    echo -e "${RED}Error:${NC} Directory '${WORK_DIR}' does not exist"
    exit 1
fi
WORK_DIR="$(realpath "${WORK_DIR}")"

HUB_REPO="bactflow"
LOCAL_IMAGE="bactflow:v1.0"
DOCKERFILE="${SCRIPT_DIR}/UI/bactflow_app/Dockerfile"
BUILD_CONTEXT="${SCRIPT_DIR}"
FALLBACK_TAG="v1.0"
DEFAULT_CPUS="10"
DEFAULT_MEMORY="16g"
CONTAINER_PORT="5000"
DEFAULT_HOST_PORT="5000"
APP_PY="${SCRIPT_DIR}/UI/bactflow_app/app.py"
APP_IN_IMAGE="/app/UI/bactflow_app/app.py"

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
    MEM_VALUE="$(normalize_memory_value "${USER_MEMORY}")"
elif [[ -n "${DEFAULT_MEMORY}" ]]; then
    MEM_VALUE="$(normalize_memory_value "${DEFAULT_MEMORY}")"
else
    MEM_VALUE=""
fi

if [[ -n "${MEM_VALUE}" ]]; then
    REQ_G="$(echo "${MEM_VALUE}" | sed -E 's/[GgMmBb]+$//')"
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

case "${SERVICE}" in
    docker) run_docker ;;
    apptainer) run_apptainer ;;
    conda) run_conda ;;
    *)
        echo -e "${RED}Error:${NC} Unknown service '${SERVICE}'"
        usage
        ;;
esac
