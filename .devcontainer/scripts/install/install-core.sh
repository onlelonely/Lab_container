#!/bin/bash

# Core packages installation script with error handling
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
readonly SCRIPT_NAME="install-core.sh"

# Environment variables with defaults
export PYTHON_VERSION="${PYTHON_VERSION:-3.11}"
export R_VERSION="${R_VERSION:-4.3}"
export WORKSPACE_DIR="${WORKSPACE_DIR:-/workspace}"
export PARALLEL_JOBS="${PARALLEL_JOBS:-4}"

# Main installation function
main() {
    log_script_start "$SCRIPT_NAME"
    
    # Parse command line arguments
    local use_parallel="${1:-false}"
    local max_jobs="${2:-$PARALLEL_JOBS}"
    
    # Wait for network connectivity
    if ! wait_for_network 60 5; then
        log_error "Network connectivity required for package installation"
        exit 1
    fi
    
    # Validate environment
    if ! validate_environment; then
        log_error "Environment validation failed"
        exit 1
    fi
    
    # Update package lists
    log_info "Updating package lists..."
    if ! retry_with_backoff 3 5 2 "apt-get update"; then
        log_error "Failed to update package lists"
        exit 1
    fi
    
    # Choose installation method
    if [ "$use_parallel" = "true" ]; then
        log_info "Using parallel installation with $max_jobs jobs"
        install_core_parallel "$max_jobs"
    else
        log_info "Using sequential installation"
        install_core_sequential
    fi
    
    # Validate installation
    log_info "Validating core installation..."
    if ! validate_python "$PYTHON_VERSION"; then
        log_warning "Python validation failed"
    fi
    
    if ! validate_r "$R_VERSION"; then
        log_warning "R validation failed"
    fi
    
    if ! validate_conda; then
        log_warning "Conda validation failed"
    fi
    
    # Test core packages
    log_info "Testing core packages..."
    test_core_packages
    
    log_script_end "$SCRIPT_NAME" 0
}

# Sequential installation (original method)
install_core_sequential() {
    # Install core conda environment
    log_info "Installing core conda environment..."
    if ! retry_with_backoff 5 10 2 "conda env update -f $CONFIG_DIR/conda/environment-core.yml"; then
        log_error "Failed to install core conda environment"
        exit 1
    fi
    
    # Install core Python packages
    log_info "Installing core Python packages..."
    if ! retry_with_backoff 5 10 2 "/opt/conda/bin/pip install -r $CONFIG_DIR/python/requirements-core.txt"; then
        log_error "Failed to install core Python packages"
        exit 1
    fi
    
    # Install core R packages
    log_info "Installing core R packages..."
    if ! retry_with_backoff 5 15 2 "/opt/conda/bin/Rscript $CONFIG_DIR/r/core-packages.R"; then
        log_error "Failed to install core R packages"
        exit 1
    fi
}

# Parallel installation (new method)
install_core_parallel() {
    local max_jobs="$1"
    
    # Optimize system for parallel operations
    optimize_for_parallel "$max_jobs"
    
    # Install core conda environment (sequential, foundation)
    log_info "Installing core conda environment..."
    if ! create_conda_env_parallel "$CONFIG_DIR/conda/environment-core.yml" "$max_jobs"; then
        log_error "Failed to install core conda environment"
        exit 1
    fi
    
    # Install Python and R packages in parallel
    log_info "Installing core Python and R packages in parallel..."
    local pids=()
    
    # Install Python packages in background
    install_pip_parallel "$CONFIG_DIR/python/requirements-core.txt" "$max_jobs" &
    pids+=($!)
    
    # Install R packages in background
    install_r_packages_parallel "$CONFIG_DIR/r/core-packages.R" "$max_jobs" &
    pids+=($!)
    
    # Wait for all background jobs
    local failed_jobs=0
    for pid in "${pids[@]}"; do
        if ! wait "$pid"; then
            failed_jobs=$((failed_jobs + 1))
        fi
    done
    
    if [ $failed_jobs -gt 0 ]; then
        log_error "$failed_jobs core package installation jobs failed"
        exit 1
    fi
}

# Test core packages function
test_core_packages() {
    local python_packages=("numpy" "pandas" "jupyter" "matplotlib")
    local r_packages=("tidyverse" "ggplot2" "data.table")
    
    # Test Python packages
    for pkg in "${python_packages[@]}"; do
        if ! validate_package "python" "$pkg"; then
            log_warning "Core Python package validation failed: $pkg"
        fi
    done
    
    # Test R packages
    for pkg in "${r_packages[@]}"; do
        if ! validate_package "r" "$pkg"; then
            log_warning "Core R package validation failed: $pkg"
        fi
    done
}

# Error handling
trap 'log_script_end "$SCRIPT_NAME" $?' EXIT

# Run main function
main "$@"