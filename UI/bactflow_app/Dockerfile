# syntax=docker/dockerfile:1
# Unified BactFlow image — multi-stage so the final layer only gets runtime bits.
#
# Stages (built in parallel where independent):
#   env-main  — conda pipeline env + Medaka/torch, then strip non-runtime files
#   env-mlst  — PubMLST sidecar, then strip duplicate BLAST/JVM/docs
#   jre       — custom jlink JRE (Nextflow / Pilon)
#   blast     — NCBI BLAST+ binaries only (not conda BLAST)
#   app       — pruned UI source (no Dockerfiles / env yml)
#   runtime   — debian-slim + copies from the stages above
#
# Build from the repository root:
#   docker build -f UI/bactflow_app/Dockerfile -t bactflow:v1.0 .
#   # or via launcher (rebuilds when needed):
#   ./bactflow.sh run /path/to/work --service docker --local --rebuild --cpus 8 --memory 20g
#
# Profiles (BACTFLOW_PROFILE):
#   full (default) — includes geNomad (TensorFlow + MMseqs2)
#   slim           — no geNomad (~half size; plasmid/virus screen unavailable)
#   docker build --build-arg BACTFLOW_PROFILE=slim -t bactflow:slim .
#
# Runtime env injected by bactflow.sh (session hard caps):
#   BACTFLOW_CPUS, BACTFLOW_MEMORY / BACTFLOW_MEMORY_GB, BACTFLOW_LAUNCH_FILE
#   Entrypoint writes ~/.bactflow/launch.json so the UI form + gauges + Nextflow
#   stay within --cpus / --memory for the whole session.


# ─── Main conda environment ───────────────────────────────────────────────────
FROM mambaorg/micromamba:1.5.10 AS env-main
USER root
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG BACTFLOW_PROFILE=full
SHELL ["/bin/bash", "-lc"]

RUN apt-get update \
    && apt-get install -y --no-install-recommends binutils \
    && rm -rf /var/lib/apt/lists/*

ENV CONDA_OVERRIDE_CUDA=""
COPY UI/bactflow_app/environment.docker.yml /tmp/environment.full.yml
COPY UI/bactflow_app/environment.docker.slim.yml /tmp/environment.slim.yml

RUN if [ "${BACTFLOW_PROFILE}" = "slim" ]; then \
        cp /tmp/environment.slim.yml /tmp/environment.yml; \
    else \
        cp /tmp/environment.full.yml /tmp/environment.yml; \
    fi \
    && rm -f /tmp/environment.full.yml /tmp/environment.slim.yml \
    && /bin/micromamba create -y -p /opt/conda/envs/bactflow -f /tmp/environment.yml \
    && /opt/conda/envs/bactflow/bin/pip install --no-cache-dir \
        torch --index-url https://download.pytorch.org/whl/cpu \
    && /opt/conda/envs/bactflow/bin/pip install --no-cache-dir medaka==2.2.1 \
    && micromamba clean -afy --all \
    && rm -rf /opt/conda/pkgs /tmp/environment.yml /root/.cache /tmp/pip-*

COPY UI/post_assem_app/docker/slim-env.sh /tmp/slim-env.sh
# Smoke-test after strip so a bad cleanup fails this stage early.
RUN bash /tmp/slim-env.sh /opt/conda/envs/bactflow \
    && rm -f /tmp/slim-env.sh \
    && /opt/conda/envs/bactflow/bin/python -c \
         "import flask, plotly, pandas, Bio, bakta, unicycler; import flye" \
    && /opt/conda/envs/bactflow/bin/seqkit version \
    && /opt/conda/envs/bactflow/bin/flye --version


# ─── MLST sidecar (independent of env-main → parallel build) ───────────────────
FROM mambaorg/micromamba:1.5.10 AS env-mlst
USER root
ARG MAMBA_DOCKERFILE_ACTIVATE=1
SHELL ["/bin/bash", "-lc"]

RUN apt-get update \
    && apt-get install -y --no-install-recommends binutils \
    && rm -rf /var/lib/apt/lists/*

COPY UI/bactflow_app/environment.mlst.yml /tmp/environment.mlst.yml
RUN /bin/micromamba create -y -p /opt/conda/envs/bactflow_mlst -f /tmp/environment.mlst.yml \
    && micromamba clean -afy --all \
    && rm -rf /opt/conda/pkgs /tmp/environment.mlst.yml

COPY UI/post_assem_app/docker/slim-mlst-env.sh /tmp/slim-mlst-env.sh
RUN bash /tmp/slim-mlst-env.sh /opt/conda/envs/bactflow_mlst \
    && rm -f /tmp/slim-mlst-env.sh \
    && test -x /opt/conda/envs/bactflow_mlst/bin/mlst \
    && test -x /opt/conda/envs/bactflow_mlst/bin/any2fasta \
    && env -u PERL5LIB -u PERL_LOCAL_LIB_ROOT -u CONDA_PREFIX \
        PATH=/opt/conda/envs/bactflow_mlst/bin:/usr/bin:/bin \
        /opt/conda/envs/bactflow_mlst/bin/mlst --version


# ─── Custom JRE for Nextflow / Pilon ───────────────────────────────────────────
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


# ─── NCBI BLAST+ only (no full toolkit tree in the final image) ────────────────
FROM debian:bookworm-slim AS blast
# libgomp1: NCBI blastn is linked against libgomp.so.1 (OpenMP).
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates curl libgomp1 \
    && rm -rf /var/lib/apt/lists/*
COPY UI/post_assem_app/docker/install-blast.sh /tmp/install-blast.sh
RUN bash /tmp/install-blast.sh 2.17.0 /opt/blast/bin \
    && rm -f /tmp/install-blast.sh \
    && /opt/blast/bin/blastn -version | grep -q '2.17.0'


# ─── Pruned application source ─────────────────────────────────────────────────
FROM debian:bookworm-slim AS app
COPY UI/bactflow_app /app/UI/bactflow_app
COPY UI/pre_assem_app /app/UI/pre_assem_app
COPY UI/assem_app /app/UI/assem_app
COPY UI/post_assem_app /app/UI/post_assem_app
RUN rm -rf \
        /app/UI/*/docker \
        /app/UI/*/Dockerfile \
        /app/UI/*/docker-compose.yml \
        /app/UI/*/environment.docker.yml \
        /app/UI/*/environment.docker.slim.yml \
        /app/UI/*/environment.mlst.yml \
        /app/UI/*/.dockerignore \
        /app/UI/*/__pycache__ \
        /app/UI/*/*/__pycache__ \
        /app/UI/*/.git \
    && find /app -type f -name '*.pyc' -delete \
    && ln -sfn /app/UI/bactflow_app/app.py /app/app.py


# ─── Final runtime image ───────────────────────────────────────────────────────
FROM debian:bookworm-slim AS runtime

ARG BACTFLOW_VERSION=v1.0
ARG BACTFLOW_PROFILE=full
LABEL org.opencontainers.image.title="BactFlow" \
      org.opencontainers.image.version="${BACTFLOW_VERSION}" \
      org.opencontainers.image.source="https://github.com/farhadm1990/bactflow" \
      org.opencontainers.image.description="Unified QC + Assembly + Post-assembly UI"

ENV BACTFLOW_VERSION=${BACTFLOW_VERSION} \
    BACTFLOW_PROFILE=${BACTFLOW_PROFILE} \
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

# Minimal OS libs only — no JDK, no micromamba, no build toolchain.
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

# Runtime payloads from earlier stages (no builders / JDKs / package caches).
COPY --from=jre /opt/jre /usr/lib/jvm/bactflow-java
COPY --from=blast /opt/blast/bin/ /usr/local/bin/
COPY --from=env-main /opt/conda/envs/bactflow /opt/conda/envs/bactflow
COPY --from=env-mlst /opt/conda/envs/bactflow_mlst /opt/conda/envs/bactflow_mlst
COPY UI/assem_app/docker/conda.sh /opt/conda/etc/profile.d/conda.sh
COPY UI/assem_app/docker/conda /opt/conda/bin/conda
COPY UI/assem_app/docker/entrypoint.sh /usr/local/bin/bactflow-entrypoint
COPY UI/assem_app/docker/pilon /usr/local/bin/pilon
COPY UI/post_assem_app/bactflow_runtime.sh /opt/conda/envs/bactflow/etc/conda/activate.d/zz-bactflow-runtime.sh
COPY --from=app /app /app

RUN mkdir -p /opt/conda/envs/bactflow/etc/conda/activate.d /opt/conda/condabin \
        /opt/conda/envs/bactflow_mlst/bin \
    && chmod 755 /opt/conda/bin/conda /usr/local/bin/bactflow-entrypoint /usr/local/bin/pilon \
        /opt/conda/envs/bactflow/etc/conda/activate.d/zz-bactflow-runtime.sh \
    && ln -sfn /opt/conda/bin/conda /opt/conda/condabin/conda \
    && ln -sfn /usr/local/bin/nextflow /opt/conda/envs/bactflow/bin/nextflow \
    && ln -sfn /usr/local/bin/blastn /opt/conda/envs/bactflow/bin/blastn \
    && ln -sfn /usr/local/bin/blastp /opt/conda/envs/bactflow/bin/blastp \
    && ln -sfn /usr/local/bin/makeblastdb /opt/conda/envs/bactflow/bin/makeblastdb \
    && ln -sfn /usr/local/bin/blastn /opt/conda/envs/bactflow_mlst/bin/blastn \
    && ln -sfn /usr/local/bin/blastp /opt/conda/envs/bactflow_mlst/bin/blastp \
    && ln -sfn /usr/local/bin/makeblastdb /opt/conda/envs/bactflow_mlst/bin/makeblastdb \
    && ln -sfn /app/UI/bactflow_app/app.py /app/app.py \
    && printf '%s\n' \
        'export CONDA_PREFIX=/opt/conda/envs/bactflow' \
        'export CONDA_DEFAULT_ENV=bactflow' \
        'export PATH="/opt/conda/envs/bactflow/bin:/usr/local/bin:/opt/conda/bin:$PATH"' \
        'export LD_LIBRARY_PATH="/opt/conda/envs/bactflow/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"' \
        '[ -f /opt/conda/etc/profile.d/conda.sh ] && . /opt/conda/etc/profile.d/conda.sh' \
        > /etc/profile.d/bactflow.sh \
    && chmod 644 /etc/profile.d/bactflow.sh \
    && echo '[ -f /etc/profile.d/bactflow.sh ] && . /etc/profile.d/bactflow.sh' >> /etc/bash.bashrc \
    && chmod +x \
        /app/UI/bactflow_app/*.py \
        /app/UI/pre_assem_app/*.sh \
        /app/UI/assem_app/*.sh \
        /app/UI/post_assem_app/*.sh \
        2>/dev/null || true \
    && python -c "import flask, plotly, pandas, Bio, psutil, pkg_resources, circlator, bakta, unicycler, reportlab" \
    && seqkit version \
    && flye --version \
    && spades.py --version \
    && unicycler --version \
    && blastn -version | grep -q '2.17.0' \
    && bakta --version \
    && command -v gtdbtk \
    && if [ "${BACTFLOW_PROFILE}" = "full" ]; then \
         command -v genomad \
         && python -c "import tensorflow, genomad; print('geNomad OK')"; \
       else \
         echo "slim profile: geNomad omitted"; \
       fi \
    && python -c "import medaka, torch; print('Medaka/torch OK', torch.__version__)" \
    && test -x /opt/conda/envs/bactflow_mlst/bin/mlst \
    && test -x /opt/conda/envs/bactflow_mlst/bin/blastn \
    && test -x /opt/conda/envs/bactflow_mlst/bin/any2fasta \
    && env -u PERL5LIB -u PERL_LOCAL_LIB_ROOT -u CONDA_PREFIX \
        PATH=/opt/conda/envs/bactflow_mlst/bin:/usr/local/bin:/usr/bin:/bin \
        /opt/conda/envs/bactflow_mlst/bin/mlst --version \
    && checkm -h >/dev/null \
    && quast.py --help >/dev/null \
    && /usr/lib/jvm/bactflow-java/bin/java -version \
    && java -jar /usr/local/share/pilon/pilon.jar --help >/dev/null \
    && pilon --help >/dev/null \
    && /usr/local/bin/nextflow -version 2>&1 | grep -qE 'version 24\.' \
    && test -x /app/UI/post_assem_app/mlst_run.sh \
    && test -f /app/UI/post_assem_app/genome_report.py \
    && test -e /app/app.py \
    && test -f /app/UI/bactflow_app/memory_cap.py \
    && test -f /app/UI/bactflow_app/static/js/chrome.js \
    && test -f /app/UI/pre_assem_app/static/js/pre-dashboard.js \
    && python -c "import sys; sys.path.insert(0,'/app/UI/bactflow_app'); from memory_cap import parse_and_cap_cpus, env_launch_cpus; assert parse_and_cap_cpus(99, host=12) >= 1" \
    && rm -rf /root/.nextflow/cache /tmp/* /root/.cache \
    && conda info --base | grep -qx /opt/conda \
    && conda env list --json | grep -q bactflow \
    && conda env list --json | grep -q bactflow_mlst

WORKDIR /app
ENV BACTFLOW_MODULE=ui
ENV BACTFLOW_HOST_PORT=5000
# Optional session caps (set by bactflow.sh --cpus / --memory at run time):
# ENV BACTFLOW_CPUS=
# ENV BACTFLOW_MEMORY=
# ENV BACTFLOW_MEMORY_GB=
# ENV BACTFLOW_LAUNCH_FILE=

EXPOSE 5000
HEALTHCHECK --interval=20s --timeout=8s --start-period=45s --retries=8 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:5000/', timeout=5)"

ENTRYPOINT ["/usr/local/bin/bactflow-entrypoint"]
CMD ["python", "/app/UI/bactflow_app/app.py"]
