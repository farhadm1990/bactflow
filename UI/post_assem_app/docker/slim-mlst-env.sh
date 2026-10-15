#!/usr/bin/env bash
# Shrink the bactflow_mlst sidecar: drop duplicate BLAST/JVM/docs.
# Runtime BLAST comes from /usr/local/bin (NCBI tarball in the main image).
set -euo pipefail

PREFIX="${1:-/opt/conda/envs/bactflow_mlst}"

if [ ! -d "${PREFIX}" ]; then
    echo "slim-mlst-env: missing prefix ${PREFIX}" >&2
    exit 1
fi

rm -rf \
    "${PREFIX}/lib/ncbi-blast+" \
    "${PREFIX}/lib/jvm" \
    "${PREFIX}/include" \
    "${PREFIX}/man" \
    "${PREFIX}/docs" \
    "${PREFIX}/doc" \
    "${PREFIX}/share/doc" \
    "${PREFIX}/share/man" \
    "${PREFIX}/share/info" \
    "${PREFIX}/share/locale" \
    "${PREFIX}/compiler_compat" \
    "${PREFIX}/x86_64-conda-linux-gnu" \
    "${PREFIX}/lib/cmake" \
    "${PREFIX}/lib/pkgconfig" \
    "${PREFIX}/share/gtk-doc" \
    "${PREFIX}/share/info" \
    "${PREFIX}/test"

# Conda BLAST / Java helpers — image already ships blastn + a custom JRE.
find "${PREFIX}/bin" -maxdepth 1 -type f \( \
    -name 'blast*' -o -name 'makeblastdb' -o -name 'makeprofiledb' \
    -o -name 'blastdb*' -o -name 'dustmasker' -o -name 'segmasker' \
    -o -name 'windowmasker' -o -name 'psiblast' -o -name 'rpsblast*' \
    -o -name 'deltablast' -o -name 'legacy_blast.pl' \
    -o -name 'update_blastdb.pl' -o -name 'cleanup-blastdb-volumes.py' \
    -o -name 'xtract*' -o -name 'transmute*' -o -name 'rchive*' \
    -o -name 'efetch' -o -name 'elink' -o -name 'einfo' -o -name 'esearch' \
    -o -name 'eutils' -o -name 'java' -o -name 'javac' -o -name 'jlink' \
    -o -name 'keytool' -o -name 'jar' \
    \) -delete 2>/dev/null || true

# Optional RNA tools pulled by bioperl stack — mlst does not need them.
find "${PREFIX}/bin" -maxdepth 1 -type f -name 'RNA*' -delete 2>/dev/null || true
rm -rf "${PREFIX}/share/ViennaRNA" "${PREFIX}/lib/libRNA.a" 2>/dev/null || true

# Drop unused Python in the sidecar if present (mlst is Perl-first).
rm -rf \
    "${PREFIX}/lib/python3.12/test" \
    "${PREFIX}/lib/python3.12/idlelib" \
    "${PREFIX}/lib/python3.12/tkinter" \
    "${PREFIX}/lib/python3.12/ensurepip" \
    "${PREFIX}/lib/python3.12/turtledemo" \
    "${PREFIX}/lib/python3.12/pydoc_data" \
    "${PREFIX}/lib/python3.12/config-"* 2>/dev/null || true

find "${PREFIX}" -type f \( -name '*.a' -o -name '*.la' -o -name '*.pyc' -o -name '*.js.map' \) -delete
find "${PREFIX}" -type d -name '__pycache__' -prune -exec rm -rf {} +

if command -v strip >/dev/null 2>&1; then
    find "${PREFIX}" -type f \( -name '*.so' -o -name '*.so.*' \) \
        ! -name 'libperl.so*' -print0 \
        | xargs -0 -r strip --strip-unneeded 2>/dev/null || true
    find "${PREFIX}/bin" -type f ! -name 'perl' ! -name 'perl5*' ! -name 'mlst' \
        ! -name 'any2fasta' -print0 \
        | xargs -0 -r strip --strip-unneeded 2>/dev/null || true
fi

echo "slim-mlst-env: $(du -sh "${PREFIX}" | awk '{print $1}')"
