#!/usr/bin/env bash
# Drop compiler/JDK/docs/static libs that conda ships but runtime does not need.
set -euo pipefail

PREFIX="${1:-/opt/conda/envs/bactflow}"

if [ ! -d "${PREFIX}" ]; then
    echo "slim-env: missing prefix ${PREFIX}" >&2
    exit 1
fi

rm -rf \
    "${PREFIX}/x86_64-conda-linux-gnu" \
    "${PREFIX}/aarch64-conda-linux-gnu" \
    "${PREFIX}/libexec/gcc" \
    "${PREFIX}/lib/gcc" \
    "${PREFIX}/lib/jvm" \
    "${PREFIX}/include" \
    "${PREFIX}/man" \
    "${PREFIX}/share/man" \
    "${PREFIX}/share/doc" \
    "${PREFIX}/share/info" \
    "${PREFIX}/share/locale" \
    "${PREFIX}/share/gir-1.0" \
    "${PREFIX}/share/X11" \
    "${PREFIX}/share/cups" \
    "${PREFIX}/share/alsa" \
    "${PREFIX}/share/aclocal" \
    "${PREFIX}/compiler_compat" \
    "${PREFIX}/lib/cmake" \
    "${PREFIX}/checkm_data"

rm -f "${PREFIX}/bin/"x86_64-conda-linux-gnu-* "${PREFIX}/bin/"aarch64-conda-linux-gnu-*

# Conda BLAST is ~550MB of NCBI C++ toolkit; Debian ncbi-blast+ replaces it.
# Keep lib/perl5: MUMmer nucmer (used by Circlator) is a Perl script and needs libperl.so.
rm -rf "${PREFIX}/lib/ncbi-blast+"
find "${PREFIX}/bin" -maxdepth 1 -type f \( \
    -name 'blast*' -o -name 'makeblastdb' -o -name 'makeprofiledb' \
    -o -name 'blastdb*' -o -name 'dustmasker' -o -name 'segmasker' \
    -o -name 'windowmasker' -o -name 'psiblast' -o -name 'rpsblast*' \
    -o -name 'deltablast' -o -name 'legacy_blast.pl' \
    -o -name 'update_blastdb.pl' -o -name 'cleanup-blastdb-volumes.py' \
    -o -name 'xtract*' -o -name 'transmute*' -o -name 'rchive*' \
    -o -name 'efetch' -o -name 'elink' -o -name 'einfo' -o -name 'esearch' \
    -o -name 'eutils' \
    \) -delete 2>/dev/null || true

# Pilon's conda recipe pulls a second OpenJDK; the image uses a wrapper + one JRE.
rm -f "${PREFIX}/bin/pilon" "${PREFIX}/bin/pilon.py"
rm -rf "${PREFIX}/share/pilon"* "${PREFIX}/opt/pilon"*

# Python bits unused at runtime
rm -rf \
    "${PREFIX}/lib/python3.11/test" \
    "${PREFIX}/lib/python3.11/idlelib" \
    "${PREFIX}/lib/python3.11/tkinter" \
    "${PREFIX}/lib/python3.11/ensurepip" \
    "${PREFIX}/lib/python3.11/turtledemo" \
    "${PREFIX}/lib/python3.11/pydoc_data" \
    "${PREFIX}/lib/python3.11/config-3.11"* \
    "${PREFIX}/lib/tcl8.6" \
    "${PREFIX}/lib/tk8.6"

find "${PREFIX}" -type f \( -name '*.a' -o -name '*.la' -o -name '*.a.*' -o -name '*.pyc' -o -name '*.js.map' \) -delete
find "${PREFIX}/lib" -maxdepth 1 -type f \( -name 'libncbi-*.a*' -o -name 'libicutest.so*' \) -delete
find "${PREFIX}" -type d -name '__pycache__' -prune -exec rm -rf {} +

# SPAdes isolate mode does not need the specialized HMM packs
rm -rf \
    "${PREFIX}/share/spades/coronaspades_hmms" \
    "${PREFIX}/share/spades/biosynthetic_spades_hmms" \
    "${PREFIX}/share/spades/rna_spades" 2>/dev/null || true

if command -v strip >/dev/null 2>&1; then
    find "${PREFIX}" -type f \( -name '*.so' -o -name '*.so.*' \) ! -name 'libperl.so*' -print0 \
        | xargs -0 -r strip --strip-unneeded 2>/dev/null || true
    find "${PREFIX}/bin" -type f ! -name 'perl' ! -name 'perl5*' -print0 \
        | xargs -0 -r strip --strip-unneeded 2>/dev/null || true
fi

# conda-meta is small and helps debugging; keep it.
rm -rf /opt/conda/pkgs /tmp/environment.yml
