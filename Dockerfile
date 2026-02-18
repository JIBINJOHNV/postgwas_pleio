# Use the platform flag to ensure compatibility with the base image
FROM --platform=linux/amd64 jibinjv/postgwas-pleio:1.0

# Switch to root to handle environment creation and permissions
USER root

# Install system-level dependencies for bioinformatics tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Set up the working directory
WORKDIR /app

# Copy the source code and configuration
COPY pyproject.toml README.md ./
COPY src/ ./src/

# 1. Create the environment with bcftools from bioconda
# 2. Install Python dependencies
# 3. Install the actual package
# 4. Fix permissions
RUN micromamba create -n postgwas -y \
    -c conda-forge \
    -c bioconda \
    python=3.9 pip bcftools && \
    micromamba run -n postgwas pip install rich rich-argparse && \
    micromamba run -n postgwas pip install --no-cache-dir . && \
    chown -R mambauser:mambauser /opt/conda && \
    micromamba clean --all --yes

# Switch back to the non-root user
USER mambauser

# This ensures that whenever you run a command, it's executed 
# inside the 'postgwas' environment automatically.
ENTRYPOINT ["micromamba", "run", "-n", "postgwas"]

# Default execution
CMD ["postgwas-pleio"]