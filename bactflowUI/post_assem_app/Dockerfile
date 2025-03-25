
FROM continuumio/miniconda3:latest AS builder


ENV CONDA_ENV=bactflow
ENV PATH="/opt/conda/envs/$CONDA_ENV/bin:$PATH"


WORKDIR /app


COPY config.yml /tmp/config.yml

# Install Mamba and create the Conda environment
RUN conda install -n base -c conda-forge mamba -y && \
    mamba env create -f /tmp/config.yml -y && \
    conda clean --all -y

RUN mamba shell init --shell bash --root-prefix=/opt/conda && \
echo "source /opt/conda/bin/activate $CONDA_ENV" >> ~/.bashrc

RUN conda init bash && \
    echo "source /opt/conda/bin/activate $CONDA_ENV" >> ~/.bashrc

# Install Nextflow
RUN wget -qO- https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow


FROM debian:bullseye-slim AS final


ENV CONDA_ENV=bactflow
ENV PATH="/opt/conda/bin:$PATH"

# Install dependencies 
RUN apt-get update && apt-get install -y \
    bash \
    wget \
    openjdk-11-jre \
    && rm -rf /var/lib/apt/lists/*

# Copy Conda installation and environment from builder stage
COPY --from=builder /opt/conda /opt/conda
COPY --from=builder /usr/local/bin/nextflow /usr/local/bin/nextflow

# Set working directory
WORKDIR /app


COPY static/ static/
COPY templates/ templates/
COPY test_out/ test_out/
COPY *.py *.sh *.r *.yml main.nf nextflow.config ./


RUN chmod +x *.sh


EXPOSE 5001



CMD ["/bin/bash", "-c", "source /opt/conda/bin/activate bactflow && python post_assembly.py"]
