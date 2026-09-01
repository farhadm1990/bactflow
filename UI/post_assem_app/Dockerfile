# syntax=docker/dockerfile:1

FROM condaforge/mambaforge:24.9.2-0 AS builder

ENV CONDA_ENV=bactflow
ENV BACTFLOW_NXF_VER=24.10.2
ENV NXF_VER=24.10.2
ENV MAMBA_NO_LOW_SPEED_LIMIT=1
ENV PATH="/opt/conda/envs/${CONDA_ENV}/bin:/opt/conda/bin:${PATH}"

WORKDIR /build

COPY config.yml bactflow_runtime.sh ./

RUN --mount=type=cache,target=/opt/conda/pkgs \
    mamba env create -f config.yml -y && \
    mamba clean -afy && \
    /opt/conda/envs/bactflow/bin/Rscript -e 'if (!requireNamespace("MatrixGenerics", quietly = TRUE)) BiocManager::install("MatrixGenerics", force = TRUE)' && \
    mkdir -p /opt/conda/envs/bactflow/etc/conda/activate.d && \
    cp bactflow_runtime.sh /opt/conda/envs/bactflow/etc/conda/activate.d/zz-bactflow-runtime.sh && \
    chmod +x /opt/conda/envs/bactflow/etc/conda/activate.d/zz-bactflow-runtime.sh && \
    /opt/conda/envs/bactflow/bin/java -version 2>&1 | grep -qE 'version "1[7-9]\.|version "2[0-4]\.' && \
    /opt/conda/envs/bactflow/bin/nextflow -version 2>&1 | grep -qE 'version 24\.' && \
    find /opt/conda/envs/bactflow -type f \( -name '*.a' -o -name '*.la' \) -delete && \
    find /opt/conda/envs/bactflow -type d -name '__pycache__' -prune -exec rm -rf {} + && \
    find /opt/conda/envs/bactflow -type f -name '*.pyc' -delete && \
    rm -rf /opt/conda/envs/bactflow/share/man /opt/conda/envs/bactflow/share/doc && \
    mamba install -y -c conda-forge conda-pack && \
    conda-pack -n bactflow -o /tmp/bactflow-env.tar.gz && \
    mamba clean -afy


FROM debian:bookworm-slim AS final

ENV CONDA_ENV=bactflow
ENV BACTFLOW_NXF_VER=24.10.2
ENV NXF_VER=24.10.2
ENV PYTHONUNBUFFERED=1
ENV PATH="/opt/conda/envs/bactflow/bin:${PATH}"
ENV BACTFLOW_JAVA_HOME=/usr/lib/jvm/bactflow-java

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    ca-certificates \
    openjdk-17-jre-headless \
    && ln -sfn /usr/lib/jvm/java-17-openjdk-amd64 /usr/lib/jvm/bactflow-java \
    && rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

COPY --from=builder /tmp/bactflow-env.tar.gz /tmp/bactflow-env.tar.gz

RUN mkdir -p /opt/conda/envs/bactflow && \
    tar -xzf /tmp/bactflow-env.tar.gz -C /opt/conda/envs/bactflow && \
    /opt/conda/envs/bactflow/conda-unpack && \
    rm /tmp/bactflow-env.tar.gz && \
    find /opt/conda/envs/bactflow -type d -name '__pycache__' -prune -exec rm -rf {} + && \
    /usr/lib/jvm/bactflow-java/bin/java -version && \
    /opt/conda/envs/bactflow/bin/nextflow -version 2>&1 | grep -qE 'version 24\.'

WORKDIR /app

COPY static/ static/
COPY templates/ templates/
COPY *.py *.sh *.r main.nf nextflow.config ./

RUN chmod +x *.sh 2>/dev/null || true

EXPOSE 5001

CMD ["python", "post_assembly.py"]
