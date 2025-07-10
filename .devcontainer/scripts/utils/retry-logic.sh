#!/bin/bash

# Retry logic with exponential backoff
set -euo pipefail

# Source logging functions
source "$(dirname "${BASH_SOURCE[0]}")/logging.sh"

# Default retry parameters
readonly DEFAULT_MAX_ATTEMPTS=3
readonly DEFAULT_DELAY=5
readonly DEFAULT_BACKOFF_FACTOR=2
readonly DEFAULT_MAX_DELAY=300

# Retry function with exponential backoff
retry_with_backoff() {
    local max_attempts="${1:-$DEFAULT_MAX_ATTEMPTS}"
    local delay="${2:-$DEFAULT_DELAY}"
    local backoff_factor="${3:-$DEFAULT_BACKOFF_FACTOR}"
    local command="$4"
    local max_delay="${5:-$DEFAULT_MAX_DELAY}"
    
    local attempt=1
    local current_delay="$delay"
    
    while [ $attempt -le $max_attempts ]; do
        log_info "Attempt $attempt of $max_attempts: $command"
        
        if eval "$command"; then
            log_success "Command succeeded on attempt $attempt"
            return 0
        fi
        
        if [ $attempt -lt $max_attempts ]; then
            log_warning "Command failed, retrying in ${current_delay}s..."
            sleep "$current_delay"
            
            # Calculate next delay with exponential backoff
            current_delay=$((current_delay * backoff_factor))
            
            # Cap delay at max_delay
            if [ $current_delay -gt $max_delay ]; then
                current_delay=$max_delay
            fi
        fi
        
        attempt=$((attempt + 1))
    done
    
    log_error "Command failed after $max_attempts attempts: $command"
    return 1
}

# Retry function for package installations
retry_package_install() {
    local package_manager="$1"
    local package_name="$2"
    local extra_args="${3:-}"
    
    local command
    case "$package_manager" in
        "conda"|"mamba"|"micromamba")
            command="$package_manager install -y $package_name $extra_args"
            ;;
        "pip")
            command="pip install $package_name $extra_args"
            ;;
        "apt")
            command="apt-get install -y $package_name $extra_args"
            ;;
        *)
            log_error "Unsupported package manager: $package_manager"
            return 1
            ;;
    esac
    
    retry_with_backoff 5 10 2 "$command" 120
}

# Retry function for network operations
retry_network_operation() {
    local operation="$1"
    retry_with_backoff 5 3 2 "$operation" 60
}

# Retry function with custom error handling
retry_with_custom_handler() {
    local max_attempts="$1"
    local delay="$2"
    local command="$3"
    local error_handler="$4"
    
    local attempt=1
    local current_delay="$delay"
    
    while [ $attempt -le $max_attempts ]; do
        log_info "Attempt $attempt of $max_attempts: $command"
        
        if eval "$command"; then
            log_success "Command succeeded on attempt $attempt"
            return 0
        fi
        
        # Run custom error handler
        if [ -n "$error_handler" ]; then
            log_info "Running custom error handler"
            eval "$error_handler"
        fi
        
        if [ $attempt -lt $max_attempts ]; then
            log_warning "Command failed, retrying in ${current_delay}s..."
            sleep "$current_delay"
            current_delay=$((current_delay * 2))
        fi
        
        attempt=$((attempt + 1))
    done
    
    log_error "Command failed after $max_attempts attempts: $command"
    return 1
}

# Function to test network connectivity
test_network_connectivity() {
    local host="${1:-8.8.8.8}"
    local timeout="${2:-5}"
    
    if ping -c 1 -W "$timeout" "$host" > /dev/null 2>&1; then
        return 0
    else
        return 1
    fi
}

# Function to wait for network connectivity
wait_for_network() {
    local max_wait="${1:-60}"
    local check_interval="${2:-5}"
    
    local elapsed=0
    
    log_info "Waiting for network connectivity..."
    
    while [ $elapsed -lt $max_wait ]; do
        if test_network_connectivity; then
            log_success "Network connectivity established"
            return 0
        fi
        
        sleep "$check_interval"
        elapsed=$((elapsed + check_interval))
        log_info "Waiting for network... (${elapsed}s elapsed)"
    done
    
    log_error "Network connectivity timeout after ${max_wait}s"
    return 1
}

# Parallel installation functions
install_packages_parallel() {
    local package_manager="$1"
    local package_list="$2"
    local max_jobs="${3:-4}"
    local extra_args="${4:-}"
    
    log_info "Installing packages in parallel with $max_jobs jobs"
    
    # Split package list into chunks
    local temp_dir=$(mktemp -d)
    local chunk_size=$(($(echo "$package_list" | wc -w) / max_jobs + 1))
    
    echo "$package_list" | tr ' ' '\n' | split -l "$chunk_size" - "$temp_dir/chunk_"
    
    # Function to install a chunk of packages
    install_chunk() {
        local chunk_file="$1"
        local packages=$(cat "$chunk_file" | tr '\n' ' ')
        
        if [ -n "$packages" ]; then
            log_info "Installing chunk: $packages"
            retry_package_install "$package_manager" "$packages" "$extra_args"
        fi
    }
    
    # Export function for parallel execution
    export -f install_chunk retry_package_install retry_with_backoff log_info log_error log_success
    
    # Install chunks in parallel
    local chunk_files=("$temp_dir"/chunk_*)
    local pids=()
    
    for chunk_file in "${chunk_files[@]}"; do
        if [ -f "$chunk_file" ]; then
            install_chunk "$chunk_file" &
            pids+=($!)
        fi
    done
    
    # Wait for all background jobs to complete
    local failed_jobs=0
    for pid in "${pids[@]}"; do
        if ! wait "$pid"; then
            failed_jobs=$((failed_jobs + 1))
        fi
    done
    
    # Cleanup
    rm -rf "$temp_dir"
    
    if [ $failed_jobs -eq 0 ]; then
        log_success "All parallel installations completed successfully"
        return 0
    else
        log_error "$failed_jobs parallel installation jobs failed"
        return 1
    fi
}

# Parallel conda environment creation
create_conda_env_parallel() {
    local env_file="$1"
    local max_jobs="${2:-4}"
    
    log_info "Creating conda environment with parallel downloads"
    
    # Set conda parallel download
    micromamba config --set channel_priority strict
    micromamba config --set solver classic
    
    # Use parallel solver if available
    if micromamba config --get solver 2>/dev/null | grep -q "libmamba"; then
        log_info "Using libmamba solver for faster dependency resolution"
    fi
    
    # Create environment with parallel downloads
    CONDA_PARALLEL_THREADS="$max_jobs" retry_with_backoff 5 15 2 "micromamba env create -f $env_file"
}

# Parallel pip installation with optimizations
install_pip_parallel() {
    local requirements_file="$1"
    local max_jobs="${2:-4}"
    
    log_info "Installing pip packages with parallel downloads"
    
    # Configure pip for parallel downloads
    pip config set global.parallel-downloads "$max_jobs"
    
    # Use pip with parallel processing
    retry_with_backoff 5 10 2 "pip install --no-cache-dir --parallel --jobs $max_jobs -r $requirements_file"
}

# Parallel R package installation
install_r_packages_parallel() {
    local packages_script="$1"
    local max_jobs="${2:-4}"
    
    log_info "Installing R packages with parallel processing"
    
    # Create parallel R installation script
    local parallel_script=$(mktemp)
    cat > "$parallel_script" << EOF
# Set parallel options
options(Ncpus = $max_jobs)
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Enable parallel installation
if (!require("parallel", quietly = TRUE)) {
    install.packages("parallel", repos = "https://cloud.r-project.org/")
}

# Source the original script with parallel support
source("$packages_script")
EOF
    
    retry_with_backoff 5 20 2 "Rscript $parallel_script"
    local result=$?
    
    rm -f "$parallel_script"
    return $result
}

# Comprehensive parallel installation orchestrator
install_all_parallel() {
    local config_dir="$1"
    local max_jobs="${2:-4}"
    local install_extended="${3:-false}"
    local install_bioconductor="${4:-true}"
    
    log_info "Starting parallel installation orchestrator"
    
    # Stage 1: Core conda environment (sequential, foundation)
    log_info "Stage 1: Installing core conda environment"
    create_conda_env_parallel "$config_dir/conda/environment-core.yml" "$max_jobs"
    
    # Stage 2: Parallel installation of core packages
    log_info "Stage 2: Installing core packages in parallel"
    {
        install_pip_parallel "$config_dir/python/requirements-core.txt" "$max_jobs" &
        install_r_packages_parallel "$config_dir/r/core-packages.R" "$max_jobs" &
        wait
    }
    
    # Stage 3: Extended packages (if requested)
    if [ "$install_extended" = "true" ]; then
        log_info "Stage 3: Installing extended packages in parallel"
        {
            retry_with_backoff 5 15 2 "micromamba env update -f $config_dir/conda/environment-extended.yml" &
            install_pip_parallel "$config_dir/python/requirements-extended.txt" "$max_jobs" &
            install_r_packages_parallel "$config_dir/r/extended-packages.R" "$max_jobs" &
            wait
        }
    fi
    
    # Stage 4: Bioinformatics packages
    log_info "Stage 4: Installing bioinformatics packages in parallel"
    {
        install_pip_parallel "$config_dir/python/requirements-bioinformatics.txt" "$max_jobs" &
        if [ "$install_bioconductor" = "true" ]; then
            install_r_packages_parallel "$config_dir/r/bioconductor-packages.R" "$max_jobs" &
        fi
        wait
    }
    
    log_success "Parallel installation orchestrator completed"
}

# Function to optimize system for parallel operations
optimize_for_parallel() {
    local max_jobs="${1:-4}"
    
    log_info "Optimizing system for parallel operations"
    
    # Set environment variables for parallel processing
    export MAKEFLAGS="-j$max_jobs"
    export CMAKE_BUILD_PARALLEL_LEVEL="$max_jobs"
    export CONDA_PARALLEL_THREADS="$max_jobs"
    
    # Configure pip for parallel downloads
    pip config set global.parallel-downloads "$max_jobs"
    
    # Set R options for parallel compilation
    export R_INSTALL_TAR="tar"
    export R_INSTALL_PARALLEL="$max_jobs"
    
    log_success "System optimized for $max_jobs parallel jobs"
}

# Export functions for use in other scripts
export -f retry_with_backoff retry_package_install retry_network_operation retry_with_custom_handler test_network_connectivity wait_for_network
export -f install_packages_parallel create_conda_env_parallel install_pip_parallel install_r_packages_parallel install_all_parallel optimize_for_parallel