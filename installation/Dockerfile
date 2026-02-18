# =====================================================================
# Base Image
# =====================================================================
FROM mambaorg/micromamba:1.4.2

USER root

# System dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    git wget ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# =====================================================================
# 1. MTAG ENV (Created from your validated YAML)
# =====================================================================
# Copy your exported environment file
COPY env_files/mtag_env.yml /tmp/mtag_env.yml

# Create the environment. 
# We don't need to specify channels here because they are in your YAML.
RUN micromamba create -y -f /tmp/mtag_env.yml && \
    micromamba clean --all --yes

# Clone the MTAG source code
WORKDIR /opt
RUN git clone https://github.com/JonJala/mtag.git && rm -rf mtag/.git

# =====================================================================
# 2. MAIN ENV: postgwas-pleio (Python 3.10)
# =====================================================================
RUN micromamba create -y -n pleio_env python=3.10 pip -c conda-forge && \
    micromamba clean --all --yes

# =====================================================================
# Install Your Package
# =====================================================================
WORKDIR /app
COPY . /app

RUN micromamba run -n pleio_env pip install --upgrade pip && \
    micromamba run -n pleio_env pip install -e .

# =====================================================================
# Final Config
# =====================================================================
# Path to your Python 2.7 executable and the MTAG script
ENV MTAG_PYTHON="/opt/conda/envs/mtag/bin/python"
ENV MTAG_SCRIPT="/opt/mtag/mtag.py"

RUN chown -R mambauser:mambauser /app /opt/mtag
USER mambauser

ENTRYPOINT ["micromamba", "run", "-n", "pleio_env"]
CMD ["python", "-m", "postgwas_pleio"]