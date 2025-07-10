#!/bin/bash

# Extended packages installation script with error handling
set -euo pipefail

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UTILS_DIR="$SCRIPT_DIR/../utils"
CONFIG_DIR="$SCRIPT_DIR/../../configs"

# Source utility functions
source "$UTILS_DIR/logging.sh"
source "$UTILS_DIR/retry-logic.sh"
source "$UTILS_DIR/validation.sh"

# Script name for logging
readonly SCRIPT_NAME="install-extended.sh"

# Environment variables with defaults
export PYTHON_VERSION="${PYTHON_VERSION:-3.11}"
export R_VERSION="${R_VERSION:-4.3}"
export WORKSPACE_DIR="${WORKSPACE_DIR:-/workspace}"
export INSTALL_BIOCONDUCTOR="${INSTALL_BIOCONDUCTOR:-true}"
export INSTALL_EXTENDED_PACKAGES="${INSTALL_EXTENDED_PACKAGES:-true}"

# Main installation function
main() {
    log_script_start "$SCRIPT_NAME"
    
    # Check if already initialized
    if [ -f ~/.devcontainer_extended_initialized ]; then
        log_info "Extended packages already initialized, skipping..."
        log_script_end "$SCRIPT_NAME" 0
        return 0
    fi
    
    # Wait for network connectivity
    if ! wait_for_network 60 5; then
        log_error "Network connectivity required for package installation"
        exit 1
    fi
    
    # Clean package caches
    log_info "Cleaning package caches..."
    if ! retry_with_backoff 3 5 2 "micromamba clean --all --yes"; then
        log_warning "Failed to clean micromamba cache"
    fi
    
    # Install extended conda environment
    if [ "$INSTALL_EXTENDED_PACKAGES" = "true" ]; then
        log_info "Installing extended conda environment..."
        if ! retry_with_backoff 5 15 2 "micromamba env update -f $CONFIG_DIR/conda/environment-extended.yml"; then
            log_error "Failed to install extended conda environment"
            exit 1
        fi
    fi
    
    # Install extended Python packages
    if [ "$INSTALL_EXTENDED_PACKAGES" = "true" ]; then
        log_info "Installing extended Python packages..."
        if ! retry_with_backoff 5 15 2 "pip install -r $CONFIG_DIR/python/requirements-extended.txt"; then
            log_error "Failed to install extended Python packages"
            exit 1
        fi
    fi
    
    # Install bioinformatics Python packages
    log_info "Installing bioinformatics Python packages..."
    if ! retry_with_backoff 5 15 2 "pip install -r $CONFIG_DIR/python/requirements-bioinformatics.txt"; then
        log_error "Failed to install bioinformatics Python packages"
        exit 1
    fi
    
    # Install Bioconductor packages
    if [ "$INSTALL_BIOCONDUCTOR" = "true" ]; then
        log_info "Installing Bioconductor packages..."
        if ! retry_with_backoff 5 20 2 "Rscript $CONFIG_DIR/r/bioconductor-packages.R"; then
            log_error "Failed to install Bioconductor packages"
            exit 1
        fi
    fi
    
    # Install extended R packages
    if [ "$INSTALL_EXTENDED_PACKAGES" = "true" ]; then
        log_info "Installing extended R packages..."
        if ! retry_with_backoff 5 20 2 "Rscript $CONFIG_DIR/r/extended-packages.R"; then
            log_error "Failed to install extended R packages"
            exit 1
        fi
    fi
    
    # Test extended packages
    log_info "Testing extended packages..."
    test_extended_packages
    
    # Mark as initialized
    touch ~/.devcontainer_extended_initialized
    log_success "Extended packages installation completed"
    
    log_script_end "$SCRIPT_NAME" 0
}

# Test extended packages function
test_extended_packages() {
    local python_packages=("biopython" "scikit-learn" "seaborn")
    local r_packages=("Biostrings" "DESeq2" "ggplot2")
    
    # Test Python packages
    for pkg in "${python_packages[@]}"; do
        if ! validate_package "python" "$pkg"; then
            log_warning "Extended Python package validation failed: $pkg"
        fi
    done
    
    # Test R packages (only if Bioconductor is enabled)
    if [ "$INSTALL_BIOCONDUCTOR" = "true" ]; then
        for pkg in "${r_packages[@]}"; do
            if ! validate_package "r" "$pkg"; then
                log_warning "Extended R package validation failed: $pkg"
            fi
        done
    fi
}

# Error handling
trap 'log_script_end "$SCRIPT_NAME" $?' EXIT

# Run main function
main "$@"