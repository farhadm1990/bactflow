#!/bin/bash
# bactflow_simple.sh - Simple BactFlow module runner with auto-browser

set -e

# Color codes
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m'

# Function to show usage
usage() {
    cat << EOF
${GREEN}BactFlow Module Runner${NC}
Run BactFlow pipeline modules with automatic browser opening

${YELLOW}USAGE:${NC}
    $0 <module> <work_dir> [options]

${YELLOW}MODULES:${NC}
    preassem    - Run pre-assembly (bactflow_preassem:v0.01) - Flask app on port 5000
    assem       - Run assembly (farhadm1990/bactflow_assem:v0.01)
    postassem   - Run post-assembly (farhadm1990/bactflow_postassem:v0.01)

${YELLOW}REQUIRED:${NC}
    work_dir    Working directory (mounted to container)

${YELLOW}OPTIONS:${NC}
    --no-browser        Don't automatically open browser
    --port PORT         Host port to map (default: 5000 for preassem, 5002 for assem and 5001 for postassem)
    --cpus N            Number of CPU cores (default: auto for preassem, 10 for assem/postassem)
    --memory SIZE       Memory limit (default: auto for preassem, 16g for assem/postassem)
    --help              Show this help message

${YELLOW}EXAMPLES:${NC}
    # Run pre-assembly with defaults
    $0 preassem /home/user/work_dir
    
    # Run pre-assembly with custom port
    $0 preassem /home/user/work_dir --port 5001
    
    # Run assembly with custom resources
    $0 assem /home/user/work_dir --cpus 20 --memory 32g
    
    # Run post-assembly without auto-browser
    $0 postassem /home/user/work_dir --no-browser
    
    # Full custom example
    $0 assem /home/user/work_dir --port 5002 --cpus 16 --memory 24g

EOF
    exit 1
}

# Function to open browser
open_browser() {
    local url=$1
    local delay=${2:-2}
    
    sleep "$delay"
    
    echo -e "${GREEN}Opening browser to: $url${NC}"
    
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        if command -v xdg-open > /dev/null; then
            xdg-open "$url" 2>/dev/null &
        elif command -v gnome-open > /dev/null; then
            gnome-open "$url" 2>/dev/null &
        elif command -v firefox > /dev/null; then
            firefox "$url" 2>/dev/null &
        else
            echo -e "${YELLOW}Could not detect browser. Please open: $url${NC}"
            return 1
        fi
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        open "$url" 2>/dev/null &
    elif [[ "$OSTYPE" == "cygwin" ]] || [[ "$OSTYPE" == "msys" ]]; then
        start "$url" 2>/dev/null &
    else
        echo -e "${YELLOW}Please open manually: $url${NC}"
        return 1
    fi
    
    return 0
}

# Function to detect available resources
detect_resources() {
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        AVAILABLE_CPU=$(nproc)
        AVAILABLE_MEM=$(free -g | awk '/^Mem:/{print $2}')g
    elif [[ "$OSTYPE" == "darwin"* ]]; then
        AVAILABLE_CPU=$(sysctl -n hw.ncpu)
        AVAILABLE_MEM=$(($(sysctl -n hw.memsize) / 1073741824))g
    else
        AVAILABLE_CPU=4
        AVAILABLE_MEM="8g"
    fi
}

# Parse arguments
OPEN_BROWSER=true
HOST_PORT=""
MODULE=""
WORK_DIR=""
USER_CPUS=""
USER_MEMORY=""

# Parse positional arguments and options
while [[ $# -gt 0 ]]; do
    case $1 in
        --no-browser)
            OPEN_BROWSER=false
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
        --help)
            usage
            ;;
        -*)
            echo -e "${RED}Error:${NC} Unknown option $1"
            usage
            ;;
        *)
            if [[ -z "$MODULE" ]]; then
                MODULE="$1"
            elif [[ -z "$WORK_DIR" ]]; then
                WORK_DIR="$1"
            else
                echo -e "${RED}Error:${NC} Too many arguments"
                usage
            fi
            shift
            ;;
    esac
done

# Check arguments
if [[ -z "$MODULE" ]] || [[ -z "$WORK_DIR" ]]; then
    echo -e "${RED}Error:${NC} Module and work directory are required"
    usage
fi

# Detect available resources
detect_resources

# Convert to absolute path
WORK_DIR=$(realpath "$WORK_DIR" 2>/dev/null || echo "$WORK_DIR")

# Validate module and set defaults
case $MODULE in
    preassem)
        IMAGE="farhadm1990/bactflow_preassem:v0.01"
        DEFAULT_CPUS=""
        DEFAULT_MEMORY=""
        CONTAINER_PORT="5000"
        DEFAULT_HOST_PORT="5000"
        ;;
    assem)
        IMAGE="farhadm1990/bactflow_assem:v0.01"
        DEFAULT_CPUS="10"
        DEFAULT_MEMORY="16g"
        CONTAINER_PORT="5002"
        DEFAULT_HOST_PORT="5002"
        ;;
    postassem)
        IMAGE="farhadm1990/bactflow_postassem:v0.01"
        DEFAULT_CPUS="10"
        DEFAULT_MEMORY="16g"
        CONTAINER_PORT="5001"
        DEFAULT_HOST_PORT="5001"
        ;;
    *)
        echo -e "${RED}Error:${NC} Unknown module '$MODULE'"
        usage
        ;;
esac

# Set CPU and memory (user specified > default > auto-detect)
if [[ -n "$USER_CPUS" ]]; then
    CPUS="--cpus=$USER_CPUS"
    CPU_VALUE="$USER_CPUS"
elif [[ -n "$DEFAULT_CPUS" ]]; then
    CPUS="--cpus=$DEFAULT_CPUS"
    CPU_VALUE="$DEFAULT_CPUS"
else
    CPUS=""
    CPU_VALUE="auto ($AVAILABLE_CPU cores)"
fi

if [[ -n "$USER_MEMORY" ]]; then
    MEMORY="--memory=$USER_MEMORY"
    MEM_VALUE="$USER_MEMORY"
elif [[ -n "$DEFAULT_MEMORY" ]]; then
    MEMORY="--memory=$DEFAULT_MEMORY"
    MEM_VALUE="$DEFAULT_MEMORY"
else
    MEMORY=""
    MEM_VALUE="auto ($AVAILABLE_MEM)"
fi

# Set host port
if [[ -z "$HOST_PORT" ]]; then
    HOST_PORT="$DEFAULT_HOST_PORT"
fi

# Validate work directory
if [ ! -d "$WORK_DIR" ]; then
    echo -e "${RED}Error:${NC} Directory '$WORK_DIR' does not exist"
    exit 1
fi


echo -e "\n${GREEN}=== Running $MODULE ===${NC}"
echo "Work directory: $WORK_DIR → $WORK_DIR"
echo "Image: $IMAGE"
echo "Resources: CPU=$CPU_VALUE, Memory=$MEM_VALUE"
echo "Port mapping: $HOST_PORT → $CONTAINER_PORT"

# Build port mapping
PORT_MAP="-p $HOST_PORT:$CONTAINER_PORT"

# Generate unique container name
CONTAINER_NAME="bactflow_${MODULE}_$$"

# Create a temporary file for logs
LOG_FILE="/tmp/bactflow_${CONTAINER_NAME}.log"
> "$LOG_FILE"  # Clear the log file

echo -e "\n${BLUE}Starting container...${NC}"
echo -e "${CYAN}Access URL: http://localhost:$HOST_PORT${NC}"
echo -e "${YELLOW}Press Ctrl+C to stop the container${NC}\n"

# Run container with port mapping and capture logs
docker run --rm \
    --name "$CONTAINER_NAME" \
    $CPUS \
    $MEMORY \
    $PORT_MAP \
    -v "$WORK_DIR:$WORK_DIR" \
    "$IMAGE" 2>&1 | tee "$LOG_FILE" &

# Get PID of background process
DOCKER_PID=$!

# Function to cleanup on exit
cleanup() {
    echo -e "\n${YELLOW}Stopping container...${NC}"
    docker stop "$CONTAINER_NAME" 2>/dev/null || true
    rm -f "$LOG_FILE"
    exit 0
}

# Set trap to catch Ctrl+C
trap cleanup INT TERM

# Wait a moment for container to start
sleep 3

# Open browser if enabled
if [[ "$OPEN_BROWSER" == true ]]; then
    # Try to open browser to the Flask app
    URL="http://localhost:$HOST_PORT"
    
    echo -e "\n${GREEN}Opening browser to $URL${NC}"
    
    # Check if the server is responding
    for i in {1..5}; do
        if curl -s -o /dev/null -w "%{http_code}" "$URL" 2>/dev/null | grep -q "200\|302\|403"; then
            open_browser "$URL" 1
            break
        elif curl -s -o /dev/null -w "%{http_code}" "http://127.0.0.1:$HOST_PORT" 2>/dev/null | grep -q "200\|302\|403"; then
            open_browser "http://127.0.0.1:$HOST_PORT" 1
            break
        elif [ $i -eq 5 ]; then
            echo -e "${YELLOW}Server not responding yet. Trying to open anyway...${NC}"
            open_browser "$URL" 1
        fi
        sleep 1
    done
fi

# Display logs in real-time with URL detection
echo -e "\n${BLUE}Container logs (real-time):${NC}\n"

# Tail the log file and detect URLs
tail -f "$LOG_FILE" | while IFS= read -r line; do
    echo "$line"
    
    # Detect Flask URLs in real-time
    if [[ "$OPEN_BROWSER" == true ]] && [[ -z "${BROWSER_OPENED:-}" ]]; then
        if [[ "$line" == *"Running on http://"* ]]; then
            # Extract URL from Flask output
            if [[ $line =~ (http://[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+:[0-9]+) ]]; then
                URL="${BASH_REMATCH[1]}"
                # Replacing container IP with localhost for browser
                URL=$(echo "$URL" | sed 's/[0-9]\+\.[0-9]\+\.[0-9]\+\.[0-9]\+/localhost/')
                echo -e "\n${GREEN}Detected Flask server: $URL${NC}"
                open_browser "$URL" 0
                BROWSER_OPENED=true
            fi
        fi
    fi
done &

# Wait for container to finish
wait $DOCKER_PID

echo -e "\n${GREEN}=== $MODULE completed successfully ===${NC}"
cleanup