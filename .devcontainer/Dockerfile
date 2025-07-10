# Multi-stage build for optimized container size and build efficiency

# Stage 1: Base system setup
FROM mcr.microsoft.com/devcontainers/base:ubuntu-22.04 AS base

# Build arguments
ARG USERNAME=devuser
ARG USER_UID=1000
ARG USER_GID=$USER_UID

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    wget \
    git \
    sudo \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user (handle existing GID/UID gracefully)
RUN if ! getent group $USER_GID > /dev/null 2>&1; then \
        groupadd --gid $USER_GID $USERNAME; \
    else \
        GROUP_NAME=$(getent group $USER_GID | cut -d: -f1); \
        if [ "$GROUP_NAME" != "$USERNAME" ]; then \
            groupadd $USERNAME; \
            USER_GID=$(getent group $USERNAME | cut -d: -f3); \
        fi; \
    fi \
    && if ! getent passwd $USER_UID > /dev/null 2>&1; then \
        useradd --uid $USER_UID --gid $USER_GID -m $USERNAME; \
    else \
        USER_NAME=$(getent passwd $USER_UID | cut -d: -f1); \
        if [ "$USER_NAME" != "$USERNAME" ]; then \
            useradd --gid $USER_GID -m $USERNAME; \
        fi; \
    fi \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# Stage 2: Development tools installation
FROM base AS dev-tools

# Copy utility scripts
COPY .devcontainer/scripts/utils/ /opt/devcontainer/utils/
RUN chmod +x /opt/devcontainer/utils/*.sh

# Set environment variables
ENV USER=$USERNAME
ENV HOME=/home/$USERNAME
ENV PATH=/opt/conda/bin:$PATH
ENV MAMBA_ROOT_PREFIX=/opt/conda

# Switch to non-root user
USER $USERNAME

# Create workspace directory
RUN sudo mkdir -p /workspace && sudo chown $USERNAME:$USERNAME /workspace

# Stage 3: Core packages (small, frequently changing)
FROM dev-tools AS core

# Copy core configuration files
COPY .devcontainer/configs/conda/environment-core.yml /tmp/environment-core.yml
COPY .devcontainer/configs/python/requirements-core.txt /tmp/requirements-core.txt
COPY .devcontainer/configs/r/core-packages.R /tmp/core-packages.R

# Install core packages with retry logic
RUN /bin/bash -c "source /opt/devcontainer/utils/logging.sh && \
    source /opt/devcontainer/utils/retry-logic.sh && \
    log_info 'Installing core packages' && \
    retry_with_backoff 5 10 2 'micromamba env create -f /tmp/environment-core.yml' && \
    retry_with_backoff 5 10 2 'pip install -r /tmp/requirements-core.txt' && \
    retry_with_backoff 5 15 2 'Rscript /tmp/core-packages.R' && \
    log_success 'Core packages installed successfully'"

# Clean up temporary files
RUN rm -rf /tmp/*.yml /tmp/*.txt /tmp/*.R

# Stage 4: Extended packages (large, less frequently changing)
FROM core AS extended

# Copy extended configuration files
COPY .devcontainer/configs/conda/environment-extended.yml /tmp/environment-extended.yml
COPY .devcontainer/configs/python/requirements-extended.txt /tmp/requirements-extended.txt
COPY .devcontainer/configs/python/requirements-bioinformatics.txt /tmp/requirements-bioinformatics.txt
COPY .devcontainer/configs/r/bioconductor-packages.R /tmp/bioconductor-packages.R
COPY .devcontainer/configs/r/extended-packages.R /tmp/extended-packages.R

# Install extended packages conditionally
ARG INSTALL_BIOCONDUCTOR=true
ARG INSTALL_EXTENDED_PACKAGES=false

RUN /bin/bash -c "source /opt/devcontainer/utils/logging.sh && \
    source /opt/devcontainer/utils/retry-logic.sh && \
    log_info 'Installing extended packages' && \
    if [ '$INSTALL_EXTENDED_PACKAGES' = 'true' ]; then \
        retry_with_backoff 5 15 2 'micromamba env update -f /tmp/environment-extended.yml' && \
        retry_with_backoff 5 15 2 'pip install -r /tmp/requirements-extended.txt'; \
    fi && \
    retry_with_backoff 5 15 2 'pip install -r /tmp/requirements-bioinformatics.txt' && \
    if [ '$INSTALL_BIOCONDUCTOR' = 'true' ]; then \
        retry_with_backoff 5 20 2 'Rscript /tmp/bioconductor-packages.R'; \
    fi && \
    if [ '$INSTALL_EXTENDED_PACKAGES' = 'true' ]; then \
        retry_with_backoff 5 20 2 'Rscript /tmp/extended-packages.R'; \
    fi && \
    log_success 'Extended packages installed successfully'"

# Clean up temporary files and caches with optimization
RUN /bin/bash -c "rm -rf /tmp/*.yml /tmp/*.txt /tmp/*.R && \
    source /opt/devcontainer/utils/optimization.sh && \
    optimize_image standard && \
    validate_optimization"

# Stage 5: Final runtime image
FROM extended AS final

# Copy installation and management scripts
COPY .devcontainer/scripts/ /opt/devcontainer/scripts/
RUN sudo chmod +x /opt/devcontainer/scripts/**/*.sh

# Copy test scripts
COPY .devcontainer/tests/ /opt/devcontainer/tests/
RUN sudo chmod +x /opt/devcontainer/tests/*.sh

# Set working directory
WORKDIR /workspace

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import sys; sys.exit(0)" && \
        R -e "q(status=0)" && \
        micromamba info > /dev/null

# Default command
CMD ["sleep", "infinity"]