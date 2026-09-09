# syntax=docker/dockerfile:1
# Slim post-assembly image: QUAST, Bakta, GTDB-Tk, CheckM, Medaka VC, R plots.
# Assemblers (Flye/SPAdes/Unicycler) live in the assembly image; databases stay on the host.
# Build from this directory:
#   cd UI/post_assem_app
#   docker build -f Dockerfile -t bactflow_postassem:v1.0 .

FROM mambaorg/micromamba:1.5.10 AS builder
USER root
ARG MAMBA_DOCKERFILE_ACTIVATE=1
SHELL ["/bin/bash", "-lc"]

RUN apt-get update \
    && apt-get install -y --no-install-recommends binutils \
    && rm -rf /var/lib/apt/lists/*

ENV CONDA_OVERRIDE_CUDA=""
COPY environment.docker.yml /tmp/environment.yml

RUN /bin/micromamba create -y -p /opt/conda/envs/bactflow -f /tmp/environment.yml \
    && /opt/conda/envs/bactflow/bin/pip install --no-cache-dir \
        torch --index-url https://download.pytorch.org/whl/cpu \
    && /opt/conda/envs/bactflow/bin/pip install --no-cache-dir medaka==2.2.1 \
    && micromamba clean -afy --all \
    && rm -f /tmp/environment.yml

COPY docker/slim-env.sh /tmp/slim-env.sh
RUN bash /tmp/slim-env.sh /opt/conda/envs/bactflow


FROM debian:bookworm-slim AS jre
RUN apt-get update \
    && apt-get install -y --no-install-recommends binutils openjdk-17-jdk-headless \
    && rm -rf /var/lib/apt/lists/* \
    && jlink \
        --add-modules java.base,java.compiler,java.datatransfer,java.desktop,java.instrument,java.logging,java.management,java.management.rmi,java.naming,java.net.http,java.prefs,java.rmi,java.scripting,java.security.jgss,java.security.sasl,java.sql,java.sql.rowset,java.xml,java.xml.crypto,jdk.crypto.cryptoki,jdk.crypto.ec,jdk.httpserver,jdk.management,jdk.unsupported,jdk.xml.dom,jdk.zipfs,jdk.localedata \
        --strip-debug \
        --no-header-files \
        --no-man-pages \
        --compress=2 \
        --output /opt/jre


FROM debian:bookworm-slim

ENV BACTFLOW_IN_DOCKER=1 \
    BACTFLOW_NXF_VER=24.10.2 \
    NXF_VER=24.10.2 \
    PYTHONUNBUFFERED=1 \
    CONDA_PREFIX=/opt/conda/envs/bactflow \
    CONDA_DEFAULT_ENV=bactflow \
    JAVA_HOME=/usr/lib/jvm/bactflow-java \
    BACTFLOW_JAVA_HOME=/usr/lib/jvm/bactflow-java \
    PATH=/opt/conda/envs/bactflow/bin:/opt/conda/bin:/usr/lib/jvm/bactflow-java/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        coreutils \
        curl \
        findutils \
        gzip \
        libgomp1 \
        libgfortran5 \
        ncbi-blast+ \
        procps \
    && rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/* \
              /usr/share/doc /usr/share/man /usr/share/info /usr/share/ncbi \
              /usr/share/perl /usr/share/perl5 /usr/lib/python3.11 \
              /usr/lib/python3 /usr/lib/valgrind \
    && rm -f /usr/bin/python3 /usr/bin/python3.11

COPY --from=jre /opt/jre /usr/lib/jvm/bactflow-java
COPY --from=builder /opt/conda/envs/bactflow /opt/conda/envs/bactflow
COPY docker/conda.sh /opt/conda/etc/profile.d/conda.sh
COPY docker/conda /opt/conda/bin/conda
COPY docker/entrypoint.sh /usr/local/bin/bactflow-entrypoint
COPY bactflow_runtime.sh /opt/conda/envs/bactflow/etc/conda/activate.d/zz-bactflow-runtime.sh

RUN mkdir -p /opt/conda/envs/bactflow/etc/conda/activate.d /opt/conda/condabin \
    && chmod 755 /opt/conda/bin/conda /usr/local/bin/bactflow-entrypoint \
        /opt/conda/envs/bactflow/etc/conda/activate.d/zz-bactflow-runtime.sh \
    && ln -sfn /opt/conda/bin/conda /opt/conda/condabin/conda \
    && if [ -x /opt/conda/envs/bactflow/bin/nextflow ]; then \
        ln -sfn /opt/conda/envs/bactflow/bin/nextflow /usr/local/bin/nextflow; \
    else \
        echo "ERROR: nextflow missing from conda env" >&2; exit 1; \
    fi \
    && printf '%s\n' \
        'export CONDA_PREFIX=/opt/conda/envs/bactflow' \
        'export CONDA_DEFAULT_ENV=bactflow' \
        'export PATH="/opt/conda/envs/bactflow/bin:/opt/conda/bin:$PATH"' \
        'export LD_LIBRARY_PATH="/opt/conda/envs/bactflow/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"' \
        '[ -f /opt/conda/etc/profile.d/conda.sh ] && . /opt/conda/etc/profile.d/conda.sh' \
        > /etc/profile.d/bactflow.sh \
    && chmod 644 /etc/profile.d/bactflow.sh \
    && echo '[ -f /etc/profile.d/bactflow.sh ] && . /etc/profile.d/bactflow.sh' >> /etc/bash.bashrc \
    && python -c "import flask, pandas, Bio, psutil, pkg_resources, circlator, bakta, torch, medaka; x=torch.zeros(4,4); x=x@x" \
    && seqkit version \
    && blastn -version \
    && diamond version \
    && bakta --version \
    && command -v gtdbtk \
    && checkm -h >/dev/null \
    && medaka --help >/dev/null \
    && quast.py --help >/dev/null \
    && Rscript -e 'library(ggplot2); library(dplyr); library(pheatmap); library(MatrixGenerics); library(Biostrings)' \
    && /usr/lib/jvm/bactflow-java/bin/java -version \
    && /usr/local/bin/nextflow -version 2>&1 | grep -qE 'version 24\.' \
    && rm -rf /root/.nextflow/cache /tmp/* \
    && conda info --base | grep -qx /opt/conda \
    && conda env list --json | grep -q bactflow

WORKDIR /app
COPY . /app/
RUN chmod +x /app/*.sh /usr/local/bin/bactflow-entrypoint /opt/conda/bin/conda 2>/dev/null || true

ENV BACTFLOW_MODULE=postassem
ENV BACTFLOW_HOST_PORT=5001

EXPOSE 5001
HEALTHCHECK --interval=20s --timeout=8s --start-period=45s --retries=8 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:5001/', timeout=5)"

ENTRYPOINT ["/usr/local/bin/bactflow-entrypoint"]
CMD ["python", "/app/post_assembly.py"]
