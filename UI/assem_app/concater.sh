#!/bin/bash
# Pool ONT barcode folders or PacBio per-sample movie/subread FASTQs.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec python3 "$SCRIPT_DIR/pool_reads.py" "$@"
