# syntax=docker/dockerfile:1
# Unified BactFlow image: one Flask hub wrapping QC, Assembly, and Post-assembly.
# Keep the runtime small: micromamba env → strip, custom JRE, NCBI BLAST tarball
# (not conda BLAST), CPU torch for Medaka. Databases stay on the host.
#
# Build from the repository root:
#   docker build -f UI/bactflow_app/Dockerfile -t bactflow:v1.0 .

FROM mambaorg/micromamba:1.5.10 AS builder
USER root
ARG MAMBA_DOCKERFILE_ACTIVATE=1
SHELL ["/bin/bash", "-lc"]

RUN apt-get update \
    && apt-get install -y --no-install-recommends binutils \
    && rm -rf /var/lib/apt/lists/*

ENV CONDA_OVERRIDE_CUDA=""
COPY UI/bactflow_app/environment.docker.yml /tmp/environment.yml
COPY UI/bactflow_app/environment.mlst.yml /tmp/environment.mlst.yml

# Main pipeline env (no mlst — perl-bioperl conflicts with infernal/samtools/spades).
RUN /bin/micromamba create -y -p /opt/conda/envs/bactflow -f /tmp/environment.yml \
    && /opt/conda/envs/bactflow/bin/pip install --no-cache-dir \
        torch --index-url https://download.pytorch.org/whl/cpu \
    && /opt/conda/envs/bactflow/bin/pip install --no-cache-dir medaka==2.2.1 \
    && micromamba clean -afy --all \
    && rm -f /tmp/environment.yml

# Sidecar env for PubMLST (Torsten Seemann mlst + any2fasta + scheme DBs).
RUN /bin/micromamba create -y -p /opt/conda/envs/bactflow_mlst -f /tmp/environment.mlst.yml \
    && micromamba clean -afy --all \
    && rm -f /tmp/environment.mlst.yml

COPY UI/post_assem_app/docker/slim-env.sh /tmp/slim-env.sh
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

ARG BACTFLOW_VERSION=v1.0
LABEL org.opencontainers.image.title="BactFlow" \
      org.opencontainers.image.version="${BACTFLOW_VERSION}" \
      org.opencontainers.image.source="https://github.com/farhadm1990/bactflow"

ENV BACTFLOW_VERSION=${BACTFLOW_VERSION} \
    BACTFLOW_IN_DOCKER=1 \
    BACTFLOW_NXF_VER=24.10.2 \
    NXF_VER=24.10.2 \
    PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    CONDA_PREFIX=/opt/conda/envs/bactflow \
    CONDA_DEFAULT_ENV=bactflow \
    JAVA_HOME=/usr/lib/jvm/bactflow-java \
    BACTFLOW_JAVA_HOME=/usr/lib/jvm/bactflow-java \
    PILON_JAR=/usr/local/share/pilon/pilon.jar \
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
        procps \
    && rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/* \
              /usr/share/doc /usr/share/man /usr/share/info /usr/share/ncbi \
              /usr/share/perl /usr/share/perl5 /usr/lib/python3.11 \
              /usr/lib/python3 /usr/lib/valgrind \
    && rm -f /usr/bin/python3 /usr/bin/python3.11 \
    && curl -fsSL -o /usr/local/bin/nextflow \
        "https://github.com/nextflow-io/nextflow/releases/download/v${NXF_VER}/nextflow" \
    && chmod 755 /usr/local/bin/nextflow \
    && mkdir -p /usr/local/share/pilon \
    && curl -fsSL -o /usr/local/share/pilon/pilon.jar \
        "https://github.com/broadinstitute/pilon/releases/download/v1.24/pilon-1.24.jar"

COPY UI/post_assem_app/docker/install-blast.sh /tmp/install-blast.sh
RUN bash /tmp/install-blast.sh 2.17.0 /usr/local/bin && rm -f /tmp/install-blast.sh

COPY --from=jre /opt/jre /usr/lib/jvm/bactflow-java
COPY --from=builder /opt/conda/envs/bactflow /opt/conda/envs/bactflow
COPY --from=builder /opt/conda/envs/bactflow_mlst /opt/conda/envs/bactflow_mlst
COPY UI/assem_app/docker/conda.sh /opt/conda/etc/profile.d/conda.sh
COPY UI/assem_app/docker/conda /opt/conda/bin/conda
COPY UI/assem_app/docker/entrypoint.sh /usr/local/bin/bactflow-entrypoint
COPY UI/assem_app/docker/pilon /usr/local/bin/pilon
COPY UI/post_assem_app/bactflow_runtime.sh /opt/conda/envs/bactflow/etc/conda/activate.d/zz-bactflow-runtime.sh

RUN mkdir -p /opt/conda/envs/bactflow/etc/conda/activate.d /opt/conda/condabin \
    && chmod 755 /opt/conda/bin/conda /usr/local/bin/bactflow-entrypoint /usr/local/bin/pilon \
        /opt/conda/envs/bactflow/etc/conda/activate.d/zz-bactflow-runtime.sh \
    && ln -sfn /opt/conda/bin/conda /opt/conda/condabin/conda \
    && ln -sfn /usr/local/bin/nextflow /opt/conda/envs/bactflow/bin/nextflow \
    && ln -sfn /usr/local/bin/blastn /opt/conda/envs/bactflow/bin/blastn \
    && ln -sfn /usr/local/bin/blastp /opt/conda/envs/bactflow/bin/blastp \
    && ln -sfn /usr/local/bin/makeblastdb /opt/conda/envs/bactflow/bin/makeblastdb \
    && printf '%s\n' \
        'export CONDA_PREFIX=/opt/conda/envs/bactflow' \
        'export CONDA_DEFAULT_ENV=bactflow' \
        'export PATH="/opt/conda/envs/bactflow/bin:/usr/local/bin:/opt/conda/bin:$PATH"' \
        'export LD_LIBRARY_PATH="/opt/conda/envs/bactflow/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"' \
        '[ -f /opt/conda/etc/profile.d/conda.sh ] && . /opt/conda/etc/profile.d/conda.sh' \
        > /etc/profile.d/bactflow.sh \
    && chmod 644 /etc/profile.d/bactflow.sh \
    && echo '[ -f /etc/profile.d/bactflow.sh ] && . /etc/profile.d/bactflow.sh' >> /etc/bash.bashrc \
    && python -c "import flask, plotly, pandas, Bio, psutil, pkg_resources, circlator, bakta, unicycler, reportlab" \
    && seqkit version \
    && flye --version \
    && spades.py --version \
    && unicycler --version \
    && blastn -version | grep -q '2.17.0' \
    && bakta --version \
    && command -v gtdbtk \
    && test -x /opt/conda/envs/bactflow_mlst/bin/mlst \
    && test -x /opt/conda/envs/bactflow_mlst/bin/any2fasta \
    && env -u PERL5LIB -u PERL_LOCAL_LIB_ROOT -u CONDA_PREFIX \
        PATH=/opt/conda/envs/bactflow_mlst/bin:/usr/bin:/bin \
        /opt/conda/envs/bactflow_mlst/bin/mlst --version \
    && checkm -h >/dev/null \
    && quast.py --help >/dev/null \
    && /usr/lib/jvm/bactflow-java/bin/java -version \
    && java -jar /usr/local/share/pilon/pilon.jar --help >/dev/null \
    && pilon --help >/dev/null \
    && /usr/local/bin/nextflow -version 2>&1 | grep -qE 'version 24\.' \
    && rm -rf /root/.nextflow/cache /tmp/* \
    && conda info --base | grep -qx /opt/conda \
    && conda env list --json | grep -q bactflow \
    && conda env list --json | grep -q bactflow_mlst

WORKDIR /app
COPY UI/bactflow_app /app/UI/bactflow_app
COPY UI/pre_assem_app /app/UI/pre_assem_app
COPY UI/assem_app /app/UI/assem_app
COPY UI/post_assem_app /app/UI/post_assem_app
RUN rm -rf \
        /app/UI/*/docker \
        /app/UI/*/Dockerfile \
        /app/UI/*/docker-compose.yml \
        /app/UI/*/environment.docker.yml \
        /app/UI/*/environment.mlst.yml \
        /app/UI/*/.dockerignore \
        /app/UI/*/__pycache__ \
        /app/UI/*/*/__pycache__ \
    && ln -sfn /app/UI/bactflow_app/app.py /app/app.py \
    && chmod +x \
        /app/UI/bactflow_app/*.py \
        /app/UI/pre_assem_app/*.sh \
        /app/UI/assem_app/*.sh \
        /app/UI/post_assem_app/*.sh \
        /usr/local/bin/bactflow-entrypoint \
        /opt/conda/bin/conda \
        /usr/local/bin/pilon \
        2>/dev/null || true \
    && test -x /app/UI/post_assem_app/mlst_run.sh \
    && test -f /app/UI/post_assem_app/genome_report.py \
    && test -f /app/UI/post_assem_app/mlst_report.py \
    && test -f /app/UI/post_assem_app/static/js/mlst-viz.js \
    && test -e /app/app.py

ENV BACTFLOW_MODULE=ui
ENV BACTFLOW_HOST_PORT=5000

EXPOSE 5000
HEALTHCHECK --interval=20s --timeout=8s --start-period=45s --retries=8 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:5000/', timeout=5)"

ENTRYPOINT ["/usr/local/bin/bactflow-entrypoint"]
CMD ["python", "/app/UI/bactflow_app/app.py"]
