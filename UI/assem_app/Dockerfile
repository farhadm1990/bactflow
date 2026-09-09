# syntax=docker/dockerfile:1
# Slim assembly image: Flask + Nextflow + Flye/SPAdes/Unicycler + QUAST/Circlator.
# BLAST and Pilon come from Debian / a single JRE (conda BLAST+OpenJDK add ~850MB).
# Taxonomy tools (GTDB-Tk, CheckM, Bakta) and Medaka stay out of this image.
# Build from this directory:
#   cd UI/assem_app
#   docker build -f Dockerfile -t bactflow_assem:v1.0 .

FROM mambaorg/micromamba:1.5.10 AS builder
USER root
ARG MAMBA_DOCKERFILE_ACTIVATE=1
SHELL ["/bin/bash", "-lc"]

RUN apt-get update \
    && apt-get install -y --no-install-recommends binutils \
    && rm -rf /var/lib/apt/lists/*

COPY environment.docker.yml /tmp/environment.yml

RUN /bin/micromamba create -y -p /opt/conda/envs/bactflow -f /tmp/environment.yml \
    && micromamba clean -afy --all \
    && rm -f /tmp/environment.yml

COPY docker/slim-env.sh /tmp/slim-env.sh
RUN bash /tmp/slim-env.sh /opt/conda/envs/bactflow


# Custom JRE is ~2x smaller than openjdk-17-jre-headless and enough for Nextflow + Pilon.
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
    PILON_JAR=/usr/local/share/pilon/pilon.jar \
    PATH=/opt/conda/envs/bactflow/bin:/opt/conda/bin:/usr/lib/jvm/bactflow-java/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        ca-certificates \
        coreutils \
        curl \
        findutils \
        gzip \
        libgomp1 \
        ncbi-blast+ \
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

COPY --from=jre /opt/jre /usr/lib/jvm/bactflow-java
COPY --from=builder /opt/conda/envs/bactflow /opt/conda/envs/bactflow
COPY docker/conda.sh /opt/conda/etc/profile.d/conda.sh
COPY docker/conda /opt/conda/bin/conda
COPY docker/entrypoint.sh /usr/local/bin/bactflow-entrypoint
COPY docker/pilon /usr/local/bin/pilon

RUN chmod 755 /opt/conda/bin/conda /usr/local/bin/bactflow-entrypoint /usr/local/bin/pilon \
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
    && python -c "import flask, pandas, Bio, psutil, pkg_resources, circlator, unicycler" \
    && seqkit version \
    && flye --version \
    && spades.py --version \
    && unicycler --version \
    && blastn -version \
    && makeblastdb -version \
    && /usr/lib/jvm/bactflow-java/bin/java -version \
    && java -jar /usr/local/share/pilon/pilon.jar --help >/dev/null \
    && pilon --help >/dev/null \
    && perl -e 'print "ok\n"' \
    && perl -c /opt/conda/envs/bactflow/bin/nucmer \
    && /usr/local/bin/nextflow -version 2>&1 | grep -qE 'version 24\.' \
    && rm -rf /root/.nextflow/cache /tmp/* \
    && conda info --base | grep -qx /opt/conda \
    && conda env list --json | grep -q bactflow

WORKDIR /app
COPY . /app/
RUN chmod +x /app/*.sh /usr/local/bin/bactflow-entrypoint /opt/conda/bin/conda /usr/local/bin/pilon 2>/dev/null || true

ENV BACTFLOW_MODULE=assem
ENV BACTFLOW_HOST_PORT=5002

EXPOSE 5002
HEALTHCHECK --interval=20s --timeout=8s --start-period=45s --retries=8 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:5002/', timeout=5)"

ENTRYPOINT ["/usr/local/bin/bactflow-entrypoint"]
CMD ["python", "/app/assembly.py"]
