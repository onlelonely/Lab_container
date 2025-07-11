#!/bin/bash

# Validation utilities for dev container
set -euo pipefail

# Source logging functions
source "$(dirname "${BASH_SOURCE[0]}")/logging.sh"

# Validate environment variables
validate_environment() {
    local required_vars=("PYTHON_VERSION" "R_VERSION" "WORKSPACE_DIR")
    local missing_vars=()
    
    log_info "Validating environment variables..."
    
    # Set default values if not already set
    export PYTHON_VERSION="${PYTHON_VERSION:-3.11}"
    export R_VERSION="${R_VERSION:-4.3}"
    export WORKSPACE_DIR="${WORKSPACE_DIR:-/workspace}"
    
    for var in "${required_vars[@]}"; do
        if [ -z "${!var:-}" ]; then
            missing_vars+=("$var")
        else
            log_info "✓ $var = ${!var}"
        fi
    done
    
    if [ ${#missing_vars[@]} -gt 0 ]; then
        log_error "Missing required environment variables: ${missing_vars[*]}"
        return 1
    fi
    
    log_success "Environment validation passed"
    return 0
}

# Validate directory exists and is writable
validate_directory() {
    local dir="$1"
    local description="${2:-directory}"
    
    if [ ! -d "$dir" ]; then
        log_error "$description does not exist: $dir"
        return 1
    fi
    
    if [ ! -w "$dir" ]; then
        log_error "$description is not writable: $dir"
        return 1
    fi
    
    log_success "✓ $description is valid: $dir"
    return 0
}

# Validate command exists
validate_command() {
    local command="$1"
    local description="${2:-$command}"
    
    if ! command -v "$command" > /dev/null 2>&1; then
        log_error "$description command not found: $command"
        return 1
    fi
    
    log_success "✓ $description is available: $(command -v "$command")"
    return 0
}

# Validate Python installation
validate_python() {
    local expected_version="${1:-3.11}"
    
    log_info "Validating Python installation..."
    
    # Try different Python commands
    local python_cmd=""
    if command -v python > /dev/null 2>&1; then
        python_cmd="python"
    elif command -v /opt/conda/bin/python > /dev/null 2>&1; then
        python_cmd="/opt/conda/bin/python"
    elif [ -f "/opt/conda/bin/python" ]; then
        python_cmd="/opt/conda/bin/python"
    fi
    
    if [ -z "$python_cmd" ]; then
        log_error "Python command not found"
        return 1
    fi
    
    log_success "✓ Python is available: $(command -v $python_cmd || echo $python_cmd)"
    
    local python_version=$($python_cmd --version 2>&1 | cut -d' ' -f2 | cut -d'.' -f1-2)
    
    if [ "$python_version" != "$expected_version" ]; then
        log_warning "Python version mismatch: expected $expected_version, got $python_version"
    else
        log_success "✓ Python version matches: $python_version"
    fi
    
    return 0
}

# Validate R installation
validate_r() {
    local expected_version="${1:-4.3}"
    
    log_info "Validating R installation..."
    
    # Try different R commands
    local r_cmd=""
    if command -v R > /dev/null 2>&1; then
        r_cmd="R"
    elif command -v /opt/conda/bin/R > /dev/null 2>&1; then
        r_cmd="/opt/conda/bin/R"
    elif [ -f "/opt/conda/bin/R" ]; then
        r_cmd="/opt/conda/bin/R"
    fi
    
    if [ -z "$r_cmd" ]; then
        log_error "R command not found"
        return 1
    fi
    
    log_success "✓ R is available: $(command -v $r_cmd || echo $r_cmd)"
    
    local r_version=$($r_cmd --version 2>/dev/null | head -n1 | grep -o '[0-9]\+\.[0-9]\+' | head -n1)
    
    if [ "$r_version" != "$expected_version" ]; then
        log_warning "R version mismatch: expected $expected_version, got $r_version"
    else
        log_success "✓ R version matches: $r_version"
    fi
    
    return 0
}

# Validate conda installation
validate_conda() {
    log_info "Validating conda installation..."
    
    # Try different conda commands
    local conda_cmd=""
    if command -v conda > /dev/null 2>&1; then
        conda_cmd="conda"
    elif command -v /opt/conda/bin/conda > /dev/null 2>&1; then
        conda_cmd="/opt/conda/bin/conda"
    elif [ -f "/opt/conda/bin/conda" ]; then
        conda_cmd="/opt/conda/bin/conda"
    fi
    
    if [ -z "$conda_cmd" ]; then
        log_error "conda command not found"
        return 1
    fi
    
    log_success "✓ conda is available: $(command -v $conda_cmd || echo $conda_cmd)"
    
    # Test conda functionality
    if $conda_cmd info > /dev/null 2>&1; then
        log_success "✓ conda is functional"
    else
        log_warning "conda info command failed (may need initialization)"
    fi
    
    return 0
}

# Validate package installation
validate_package() {
    local package_type="$1"
    local package_name="$2"
    
    case "$package_type" in
        "python")
            if python -c "import $package_name" 2>/dev/null; then
                log_success "✓ Python package available: $package_name"
                return 0
            else
                log_error "Python package not found: $package_name"
                return 1
            fi
            ;;
        "r")
            if R -e "library($package_name)" 2>/dev/null; then
                log_success "✓ R package available: $package_name"
                return 0
            else
                log_error "R package not found: $package_name"
                return 1
            fi
            ;;
        "conda")
            if conda list | grep -q "^$package_name"; then
                log_success "✓ Conda package available: $package_name"
                return 0
            else
                log_error "Conda package not found: $package_name"
                return 1
            fi
            ;;
        *)
            log_error "Unsupported package type: $package_type"
            return 1
            ;;
    esac
}

# Validate system resources
validate_system_resources() {
    log_info "Validating system resources..."
    
    # Check available memory
    local memory_gb=$(free -g | awk '/^Mem:/{print $2}')
    if [ "$memory_gb" -lt 4 ]; then
        log_warning "Low memory detected: ${memory_gb}GB (recommended: 4GB+)"
    else
        log_success "✓ Memory: ${memory_gb}GB"
    fi
    
    # Check available disk space
    local disk_gb=$(df -BG /workspace | awk 'NR==2 {print $4}' | sed 's/G//')
    if [ "$disk_gb" -lt 10 ]; then
        log_warning "Low disk space detected: ${disk_gb}GB (recommended: 10GB+)"
    else
        log_success "✓ Disk space: ${disk_gb}GB"
    fi
    
    return 0
}

# Comprehensive validation
validate_all() {
    local validation_failed=0
    
    log_info "Running comprehensive validation..."
    
    # Validate environment
    if ! validate_environment; then
        validation_failed=1
    fi
    
    # Validate directories
    if ! validate_directory "/workspace" "Workspace directory"; then
        validation_failed=1
    fi
    
    # Validate tools
    if ! validate_python; then
        validation_failed=1
    fi
    
    if ! validate_r; then
        validation_failed=1
    fi
    
    if ! validate_conda; then
        validation_failed=1
    fi
    
    # Validate system resources
    validate_system_resources
    
    if [ $validation_failed -eq 0 ]; then
        log_success "All validations passed"
        return 0
    else
        log_error "Some validations failed"
        return 1
    fi
}

# Export functions for use in other scripts
export -f validate_environment validate_directory validate_command validate_python validate_r validate_conda validate_package validate_system_resources validate_all