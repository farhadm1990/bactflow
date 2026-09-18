#!/usr/bin/env bash
# Bakta 1.12 requires BLAST+ 2.17.0. Debian bookworm ships 2.12; conda BLAST is huge.
# Install the official NCBI 2.17.0 binaries only.
set -euo pipefail

BLAST_VER="${1:-2.17.0}"
DEST="${2:-/usr/local/bin}"
url="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VER}/ncbi-blast-${BLAST_VER}+-x64-linux.tar.gz"
tmp="$(mktemp -d)"
trap 'rm -rf "${tmp}"' EXIT

echo "Installing NCBI BLAST+ ${BLAST_VER} from NCBI FTP..."
curl -fsSL --retry 8 --retry-delay 4 --retry-all-errors --connect-timeout 20 --max-time 180 \
    -o "${tmp}/blast.tgz" "${url}"
tar -xzf "${tmp}/blast.tgz" -C "${tmp}"
src="$(find "${tmp}" -type d -name 'ncbi-blast-*' | head -1)"
mkdir -p "${DEST}"
for bin in blastn blastp blastx makeblastdb blastdbcmd
do
    if [ -x "${src}/bin/${bin}" ]; then
        cp -f "${src}/bin/${bin}" "${DEST}/${bin}"
        chmod 755 "${DEST}/${bin}"
    fi
done
"${DEST}/blastn" -version | grep -q "${BLAST_VER}"
echo "blastn $("${DEST}/blastn" -version 2>&1 | head -1)"
