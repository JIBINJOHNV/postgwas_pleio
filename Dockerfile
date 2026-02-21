############################################################
# ===================== BUILDER STAGE =====================
############################################################
FROM --platform=linux/amd64 jibinjv/postgwas-pleio:1.0 AS builder

USER root

# ---------------------------------------------------------
# System deps for compilation ONLY
# ---------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps git build-essential gfortran \
    libcurl4-openssl-dev libssl-dev libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# ---------------------------------------------------------
# Clone required tools
# ---------------------------------------------------------
RUN git clone --depth 1 https://github.com/cuelee/pleio.git /opt/pleio && \
    git clone --depth 1 https://github.com/bulik/ldsc.git /opt/pleio/ldsc

# ---------------------------------------------------------
# Create environments (FULL BUILD)
# ---------------------------------------------------------
RUN micromamba create -n postgwas -y -c conda-forge -c bioconda \
    python=3.11 pip bcftools \
    polars=0.20.21 pyarrow=14 \
    pandas numpy scipy bitarray pyyaml poetry \
    compilers make \
    r-base=4.3.3 r-remotes r-data.table r-matrix r-tidyverse \
    r-future r-parallelly r-rcpp r-rcpparmadillo r-plyr r-argparse \
    r-reticulate r-rcpptoml \
    r-biocmanager r-lavaan r-psych r-igraph r-qgraph && \
    micromamba create -n pleio -y -c conda-forge \
    python=3.8 pandas=1.1.5 numpy=1.21 scipy=1.7 bitarray pyyaml

# ---------------------------------------------------------
# Copy project
# ---------------------------------------------------------
COPY pyproject.toml README.md ./
COPY src/ ./src/

# ---------------------------------------------------------
# Install Python + heavy R packages
# ---------------------------------------------------------
RUN micromamba run -n postgwas pip install rich rich-argparse && \
    micromamba run -n postgwas pip install --no-deps --no-cache-dir . && \
    micromamba run -n postgwas Rscript -e "options(Ncpus=4); remotes::install_github('sbstatgen/ASSET', upgrade='never')" && \
    micromamba run -n postgwas Rscript -e "options(Ncpus=4); remotes::install_github('gqi/fastASSET', upgrade='never')" && \
    micromamba run -n postgwas Rscript -e "options(Ncpus=4); remotes::install_github('GenomicSEM/GenomicSEM', upgrade='never')"

# ---------------------------------------------------------
# CLEAN BUILDER ENV (important before copy)
# ---------------------------------------------------------
RUN rm -rf /opt/pleio/.git /opt/pleio/ldsc/.git && \
    micromamba clean --all --yes && \
    rm -rf /opt/conda/pkgs && \
    find /opt/conda -name "*.a" -delete && \
    find /opt/conda -name "*.h" -delete && \
    find /opt/conda -type d -name "tests" -exec rm -rf {} + && \
    find /opt/conda -type d -name "__pycache__" -exec rm -rf {} +


############################################################
# ===================== RUNTIME STAGE =====================
############################################################
FROM --platform=linux/amd64 mambaorg/micromamba:1.5.8

USER root

# ---------------------------------------------------------
# Runtime minimal deps only
# ---------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps libcurl4-openssl-dev libssl-dev libxml2 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# ---------------------------------------------------------
# Copy ONLY built environments + tools
# ---------------------------------------------------------
COPY --from=builder /opt/conda /opt/conda
COPY --from=builder /opt/pleio /opt/pleio
COPY --from=builder /app /app

# ---------------------------------------------------------
# Final cleanup + ownership
# ---------------------------------------------------------
RUN find /opt/conda -type f -executable -exec strip --strip-unneeded {} + || true && \
    chown -R mambauser:mambauser /opt/conda /opt/pleio /app

ENV PATH="/opt/pleio:/opt/pleio/ldsc:${PATH}"

USER mambauser

ENTRYPOINT ["micromamba", "run", "-n", "postgwas"]
CMD ["postgwas-pleio"]