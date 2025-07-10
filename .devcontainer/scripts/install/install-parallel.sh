#!/bin/bash

# Parallel installation orchestrator script
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
readonly SCRIPT_NAME="install-parallel.sh"

# Environment variables with defaults
export PYTHON_VERSION="${PYTHON_VERSION:-3.11}"
export R_VERSION="${R_VERSION:-4.3}"
export WORKSPACE_DIR="${WORKSPACE_DIR:-/workspace}"
export INSTALL_BIOCONDUCTOR="${INSTALL_BIOCONDUCTOR:-true}"
export INSTALL_EXTENDED_PACKAGES="${INSTALL_EXTENDED_PACKAGES:-false}"
export PARALLEL_JOBS="${PARALLEL_JOBS:-4}"

# Main installation function
main() {
    log_script_start "$SCRIPT_NAME"
    
    # Parse command line arguments
    local install_type="${1:-all}"
    local max_jobs="${2:-$PARALLEL_JOBS}"
    
    # Validate environment
    if ! validate_environment; then
        log_error "Environment validation failed"
        exit 1
    fi
    
    # Wait for network connectivity
    if ! wait_for_network 60 5; then
        log_error "Network connectivity required for package installation"
        exit 1
    fi
    
    # Optimize system for parallel operations
    optimize_for_parallel "$max_jobs"
    
    # Determine CPU cores and adjust parallel jobs if needed
    local cpu_cores=$(nproc)
    if [ "$max_jobs" -gt "$cpu_cores" ]; then
        log_warning "Reducing parallel jobs from $max_jobs to $cpu_cores (available CPU cores)"
        max_jobs=$cpu_cores
    fi
    
    log_info "Starting parallel installation with $max_jobs jobs"
    
    case "$install_type" in
        "core")
            install_core_parallel "$max_jobs"
            ;;
        "extended")
            install_extended_parallel "$max_jobs"
            ;;
        "bioinformatics")
            install_bioinformatics_parallel "$max_jobs"
            ;;
        "all")
            install_all_parallel "$CONFIG_DIR" "$max_jobs" "$INSTALL_EXTENDED_PACKAGES" "$INSTALL_BIOCONDUCTOR"
            ;;
        *)
            log_error "Unknown installation type: $install_type"
            log_info "Available types: core, extended, bioinformatics, all"
            exit 1
            ;;
    esac
    
    # Validate installation
    log_info "Validating parallel installation..."
    validate_parallel_installation "$install_type"
    
    log_script_end "$SCRIPT_NAME" 0
}

# Core packages parallel installation
install_core_parallel() {
    local max_jobs="$1"
    
    log_info "Installing core packages in parallel"
    
    # Stage 1: Core conda environment
    create_conda_env_parallel "$CONFIG_DIR/conda/environment-core.yml" "$max_jobs"
    
    # Stage 2: Parallel installation of Python and R core packages
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
    
    if [ $failed_jobs -eq 0 ]; then
        log_success "Core packages installation completed successfully"
    else
        log_error "$failed_jobs core package installation jobs failed"
        exit 1
    fi
}

# Extended packages parallel installation
install_extended_parallel() {
    local max_jobs="$1"
    
    log_info "Installing extended packages in parallel"
    
    local pids=()
    
    # Update conda environment
    retry_with_backoff 5 15 2 "micromamba env update -f $CONFIG_DIR/conda/environment-extended.yml" &
    pids+=($!)
    
    # Install extended Python packages in background
    install_pip_parallel "$CONFIG_DIR/python/requirements-extended.txt" "$max_jobs" &
    pids+=($!)
    
    # Install extended R packages in background
    install_r_packages_parallel "$CONFIG_DIR/r/extended-packages.R" "$max_jobs" &
    pids+=($!)
    
    # Wait for all background jobs
    local failed_jobs=0
    for pid in "${pids[@]}"; do
        if ! wait "$pid"; then
            failed_jobs=$((failed_jobs + 1))
        fi
    done
    
    if [ $failed_jobs -eq 0 ]; then
        log_success "Extended packages installation completed successfully"
    else
        log_error "$failed_jobs extended package installation jobs failed"
        exit 1
    fi
}

# Bioinformatics packages parallel installation
install_bioinformatics_parallel() {
    local max_jobs="$1"
    
    log_info "Installing bioinformatics packages in parallel"
    
    local pids=()
    
    # Install bioinformatics Python packages in background
    install_pip_parallel "$CONFIG_DIR/python/requirements-bioinformatics.txt" "$max_jobs" &
    pids+=($!)
    
    # Install Bioconductor packages in background (if enabled)
    if [ "$INSTALL_BIOCONDUCTOR" = "true" ]; then
        install_r_packages_parallel "$CONFIG_DIR/r/bioconductor-packages.R" "$max_jobs" &
        pids+=($!)
    fi
    
    # Wait for all background jobs
    local failed_jobs=0
    for pid in "${pids[@]}"; do
        if ! wait "$pid"; then
            failed_jobs=$((failed_jobs + 1))
        fi
    done
    
    if [ $failed_jobs -eq 0 ]; then
        log_success "Bioinformatics packages installation completed successfully"
    else
        log_error "$failed_jobs bioinformatics package installation jobs failed"
        exit 1
    fi
}

# Validate parallel installation
validate_parallel_installation() {
    local install_type="$1"
    
    log_info "Validating $install_type installation"
    
    case "$install_type" in
        "core"|"all")
            # Test core packages
            local core_python_packages=("numpy" "pandas" "jupyter" "matplotlib")
            local core_r_packages=("tidyverse" "ggplot2" "data.table")
            
            for pkg in "${core_python_packages[@]}"; do
                if ! validate_package "python" "$pkg"; then
                    log_warning "Core Python package validation failed: $pkg"
                fi
            done
            
            for pkg in "${core_r_packages[@]}"; do
                if ! validate_package "r" "$pkg"; then
                    log_warning "Core R package validation failed: $pkg"
                fi
            done
            ;;
    esac
    
    case "$install_type" in
        "extended"|"all")
            # Test extended packages
            if [ "$INSTALL_EXTENDED_PACKAGES" = "true" ]; then
                local extended_python_packages=("scikit-learn" "seaborn" "plotly")
                local extended_r_packages=("shiny" "plotly" "DT")
                
                for pkg in "${extended_python_packages[@]}"; do
                    if ! validate_package "python" "$pkg"; then
                        log_warning "Extended Python package validation failed: $pkg"
                    fi
                done
                
                for pkg in "${extended_r_packages[@]}"; do
                    if ! validate_package "r" "$pkg"; then
                        log_warning "Extended R package validation failed: $pkg"
                    fi
                done
            fi
            ;;
    esac
    
    case "$install_type" in
        "bioinformatics"|"all")
            # Test bioinformatics packages
            local bio_python_packages=("biopython" "pysam")
            local bio_r_packages=("Biostrings" "GenomicRanges")
            
            for pkg in "${bio_python_packages[@]}"; do
                if ! validate_package "python" "$pkg"; then
                    log_warning "Bioinformatics Python package validation failed: $pkg"
                fi
            done
            
            if [ "$INSTALL_BIOCONDUCTOR" = "true" ]; then
                for pkg in "${bio_r_packages[@]}"; do
                    if ! validate_package "r" "$pkg"; then
                        log_warning "Bioinformatics R package validation failed: $pkg"
                    fi
                done
            fi
            ;;
    esac
    
    log_success "Parallel installation validation completed"
}

# Show usage information
usage() {
    cat << EOF
Usage: $0 [TYPE] [JOBS]

Install packages using parallel processing.

Arguments:
  TYPE    Installation type (default: all)
          - core: Core packages only
          - extended: Extended packages only
          - bioinformatics: Bioinformatics packages only
          - all: All packages

  JOBS    Number of parallel jobs (default: $PARALLEL_JOBS)
          Will be automatically adjusted based on available CPU cores

Environment Variables:
  INSTALL_BIOCONDUCTOR=true|false      Install Bioconductor packages
  INSTALL_EXTENDED_PACKAGES=true|false Install extended packages
  PARALLEL_JOBS=N                      Default number of parallel jobs

Examples:
  $0                          # Install all packages with default parallelism
  $0 core                     # Install core packages only
  $0 extended 8               # Install extended packages with 8 parallel jobs
  $0 bioinformatics 2         # Install bioinformatics packages with 2 parallel jobs

EOF
}

# Check for help flag
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

# Error handling
trap 'log_script_end "$SCRIPT_NAME" $?' EXIT

# Run main function
main "$@"