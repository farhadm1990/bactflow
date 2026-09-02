# syntax=docker/dockerfile:1
# Slim assembly image: Flask + Nextflow + Flye/SPAdes/Unicycler + QUAST/Circlator.
# Taxonomy tools (GTDB-Tk, CheckM, Bakta) and Medaka stay out of this image.
# Build from this directory:
#   cd UI/assem_app
#   docker build -f Dockerfile -t bactflow_assem:v1.0 .

FROM mambaorg/micromamba:1.5.10 AS builder
USER root
ARG MAMBA_DOCKERFILE_ACTIVATE=1
SHELL ["/bin/bash", "-lc"]

COPY environment.docker.yml /tmp/environment.yml

RUN /bin/micromamba create -y -p /opt/conda/envs/bactflow -f /tmp/environment.yml \
    && micromamba clean -afy --all \
    && find /opt/conda/envs/bactflow -type f \( -name '*.a' -o -name '*.la' \) -delete \
    && find /opt/conda/envs/bactflow -type d -name '__pycache__' -prune -exec rm -rf {} + \
    && find /opt/conda/envs/bactflow -type f -name '*.pyc' -delete \
    && rm -rf /opt/conda/envs/bactflow/share/man \
              /opt/conda/envs/bactflow/share/doc \
              /opt/conda/pkgs \
              /tmp/environment.yml


FROM debian:bookworm-slim

ENV BACTFLOW_IN_DOCKER=1 \
    BACTFLOW_NO_BROWSER=1 \
    BACTFLOW_NXF_VER=24.10.2 \
    NXF_VER=24.10.2 \
    PYTHONUNBUFFERED=1 \
    CONDA_PREFIX=/opt/conda/envs/bactflow \
    CONDA_DEFAULT_ENV=bactflow \
    JAVA_HOME=/usr/lib/jvm/bactflow-java \
    BACTFLOW_JAVA_HOME=/usr/lib/jvm/bactflow-java \
    PATH=/opt/conda/envs/bactflow/bin:/opt/conda/bin:/usr/lib/jvm/bactflow-java/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        coreutils \
        curl \
        findutils \
        git \
        gzip \
        libgomp1 \
        openjdk-17-jre-headless \
        procps \
    && ln -sfn /usr/lib/jvm/java-17-openjdk-amd64 /usr/lib/jvm/bactflow-java \
    && curl -fsSL "https://github.com/nextflow-io/nextflow/releases/download/v${NXF_VER}/nextflow" \
        -o /usr/local/bin/nextflow \
    && chmod 755 /usr/local/bin/nextflow \
    && /usr/local/bin/nextflow -version \
    && rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

COPY --from=builder /opt/conda/envs/bactflow /opt/conda/envs/bactflow
COPY docker/conda.sh /opt/conda/etc/profile.d/conda.sh
COPY docker/conda /opt/conda/bin/conda
COPY docker/entrypoint.sh /usr/local/bin/bactflow-entrypoint

RUN chmod 755 /opt/conda/bin/conda /usr/local/bin/bactflow-entrypoint \
    && mkdir -p /opt/conda/condabin \
    && ln -sfn /opt/conda/bin/conda /opt/conda/condabin/conda \
    && printf '%s\n' \
        'export CONDA_PREFIX=/opt/conda/envs/bactflow' \
        'export CONDA_DEFAULT_ENV=bactflow' \
        'export PATH="/opt/conda/envs/bactflow/bin:/opt/conda/bin:$PATH"' \
        'export LD_LIBRARY_PATH="/opt/conda/envs/bactflow/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"' \
        '[ -f /opt/conda/etc/profile.d/conda.sh ] && . /opt/conda/etc/profile.d/conda.sh' \
        > /etc/profile.d/bactflow.sh \
    && chmod 644 /etc/profile.d/bactflow.sh \
    && echo '[ -f /etc/profile.d/bactflow.sh ] && . /etc/profile.d/bactflow.sh' >> /etc/bash.bashrc \
    && python -c "import flask, pandas, Bio, psutil, pkg_resources, circlator" \
    && seqkit version \
    && flye --version \
    && spades.py --version \
    && unicycler --version \
    && /usr/lib/jvm/bactflow-java/bin/java -version \
    && nextflow -version 2>&1 | grep -qE 'version 24\.' \
    && conda info --base | grep -qx /opt/conda \
    && conda env list --json | grep -q bactflow

WORKDIR /app
COPY . /app/
RUN chmod +x /app/*.sh /usr/local/bin/bactflow-entrypoint /opt/conda/bin/conda 2>/dev/null || true

EXPOSE 5002
HEALTHCHECK --interval=20s --timeout=8s --start-period=45s --retries=8 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:5002/', timeout=5)"

ENTRYPOINT ["/usr/local/bin/bactflow-entrypoint"]
CMD ["python", "/app/assembly.py"]
