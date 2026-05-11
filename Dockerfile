############################################################
# ===================== BUILDER STAGE =====================
############################################################

# Default to linux/amd64 for Google Batch / Nextflow cloud execution.
# You can still override this using:
# docker buildx build --platform linux/amd64 ...
ARG TARGETPLATFORM=linux/amd64

FROM --platform=$TARGETPLATFORM jibinjv/postgwas-pleio:1.0 AS builder

USER root

# ---------------------------------------------------------
# System deps for compilation/building ONLY
# ---------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    git \
    wget \
    ca-certificates \
    build-essential \
    gfortran \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# ---------------------------------------------------------
# MTAG ENV
# Correct YAML location:
# installation/env_files/mtag_env.yml
# ---------------------------------------------------------
COPY installation/env_files/mtag_env.yml /tmp/mtag_env.yml

RUN micromamba create -y -f /tmp/mtag_env.yml && \
    micromamba clean --all --yes && \
    rm -f /tmp/mtag_env.yml

# ---------------------------------------------------------
# Clone required tools
# Remove existing folders first because the base image may
# already contain /opt/mtag or other tool directories.
# ---------------------------------------------------------
RUN rm -rf /opt/pleio /opt/placo /opt/mtag && \
    git clone --depth 1 https://github.com/cuelee/pleio.git /opt/pleio && \
    git clone --depth 1 https://github.com/bulik/ldsc.git /opt/pleio/ldsc && \
    git clone --depth 1 https://github.com/RayDebashree/PLACO.git /opt/placo && \
    git clone --depth 1 https://github.com/JonJala/mtag.git /opt/mtag

RUN rm -f /opt/placo/PLACO_v0.2.0_example.R && \
    rm -f /opt/placo/PLACO_v0.2.0_manual.pdf && \
    rm -rf /opt/pleio/ldsc/test

# ---------------------------------------------------------
# Create environments
# MTAG env is already created above from the YAML file.
# ---------------------------------------------------------
RUN micromamba create -n postgwas -y \
    -c conda-forge -c bioconda -c dnachun \
    python=3.11 \
    pip \
    bcftools \
    polars=0.20.21 \
    pyarrow=14 \
    psutil \
    pandas \
    numpy \
    scipy \
    bitarray \
    pyyaml \
    poetry \
    metal \
    compilers \
    make \
    r-base=4.3.3 \
    r-remotes \
    r-data.table \
    r-matrix \
    r-tidyverse \
    r-future \
    r-parallelly \
    r-rcpp \
    r-rcpparmadillo \
    r-plyr \
    r-argparse \
    r-reticulate \
    r-rcpptoml \
    r-checkmate \
    r-biocmanager \
    r-lavaan \
    r-psych \
    r-igraph \
    r-qgraph && \
    micromamba create -n pleio -y \
    -c conda-forge \
    python=3.8 \
    pandas=1.1.5 \
    numpy=1.21 \
    scipy=1.7 \
    bitarray \
    pyyaml && \
    micromamba clean --all --yes

# ---------------------------------------------------------
# Copy project
# ---------------------------------------------------------
COPY pyproject.toml README.md ./
COPY src/ ./src/

# ---------------------------------------------------------
# Install Python package + heavy R packages
# ---------------------------------------------------------
RUN micromamba run -n postgwas pip install --upgrade pip && \
    micromamba run -n postgwas pip install rich rich-argparse && \
    micromamba run -n postgwas pip install --no-deps --no-cache-dir . && \
    micromamba run -n postgwas Rscript -e "options(Ncpus=4); remotes::install_github('sbstatgen/ASSET', upgrade='never')" && \
    micromamba run -n postgwas Rscript -e "options(Ncpus=4); remotes::install_github('gqi/fastASSET', upgrade='never')" && \
    micromamba run -n postgwas Rscript -e "options(Ncpus=4); remotes::install_github('GenomicSEM/GenomicSEM', upgrade='never')"

# ---------------------------------------------------------
# Optional sanity checks
# These confirm that the expected executables/scripts exist.
# ---------------------------------------------------------
RUN test -x /opt/conda/envs/postgwas/bin/python && \
    test -x /opt/conda/envs/mtag/bin/python && \
    test -f /opt/mtag/mtag.py && \
    test -d /opt/pleio && \
    test -d /opt/pleio/ldsc && \
    test -d /opt/placo

# ---------------------------------------------------------
# Clean builder env before copying to runtime
# ---------------------------------------------------------
RUN rm -rf \
    /opt/pleio/.git \
    /opt/pleio/ldsc/.git \
    /opt/placo/.git \
    /opt/mtag/.git && \
    micromamba clean --all --yes && \
    rm -rf /opt/conda/pkgs && \
    find /opt/conda -name "*.a" -delete && \
    find /opt/conda -name "*.h" -delete && \
    find /opt/conda -type d -name "tests" -exec rm -rf {} + && \
    find /opt/conda -type d -name "__pycache__" -exec rm -rf {} +


############################################################
# ===================== RUNTIME STAGE =====================
############################################################

ARG TARGETPLATFORM=linux/amd64

FROM --platform=$TARGETPLATFORM mambaorg/micromamba:1.5.8

USER root

# ---------------------------------------------------------
# Runtime deps + Docker CLI
# Docker CLI is useful if your pipeline needs Docker-in-Docker
# or calls docker from inside the container.
# ---------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2 \
    ca-certificates \
    curl \
    gnupg \
    && install -m 0755 -d /etc/apt/keyrings \
    && curl -fsSL https://download.docker.com/linux/debian/gpg \
    | gpg --dearmor -o /etc/apt/keyrings/docker.gpg \
    && chmod a+r /etc/apt/keyrings/docker.gpg \
    && echo "deb [arch=amd64 signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/debian bullseye stable" \
    | tee /etc/apt/sources.list.d/docker.list > /dev/null \
    && apt-get update && apt-get install -y --no-install-recommends docker-ce-cli \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# ---------------------------------------------------------
# Copy built environments + tools from builder
# ---------------------------------------------------------
COPY --from=builder /opt/conda /opt/conda
COPY --from=builder /opt/pleio /opt/pleio
COPY --from=builder /opt/mtag /opt/mtag
COPY --from=builder /opt/placo /opt/placo
COPY --from=builder /app /app

# ---------------------------------------------------------
# Environment variables
# Make postgwas directly available without micromamba run.
# This is better for Nextflow / Google Batch because they
# inject the command directly.
# ---------------------------------------------------------
ENV PATH="/opt/conda/envs/postgwas/bin:/opt/pleio:/opt/pleio/ldsc:${PATH}"
ENV CONDA_PREFIX="/opt/conda/envs/postgwas"

# MTAG settings
# This assumes installation/env_files/mtag_env.yml creates env name: mtag
ENV MTAG_PYTHON="/opt/conda/envs/mtag/bin/python"
ENV MTAG_SCRIPT="/opt/mtag/mtag.py"

# Optional PLACO/Pleio helper paths
ENV PLEIO_PATH="/opt/pleio"
ENV LDSC_PATH="/opt/pleio/ldsc"
ENV PLACO_PATH="/opt/placo"

# ---------------------------------------------------------
# Final cleanup
# Keep USER root for Nextflow + Google Batch.
# ---------------------------------------------------------
RUN find /opt/conda -type f -executable -exec strip --strip-unneeded {} + || true && \
    chmod -R a+rX /opt/conda /opt/pleio /opt/mtag /opt/placo /app

USER root

ENTRYPOINT []