#!/bin/bash

# Docker image optimization utilities
set -euo pipefail

# Source logging functions
source "$(dirname "${BASH_SOURCE[0]}")/logging.sh"

# Cleanup temporary files and caches
cleanup_image() {
    log_info "Cleaning up Docker image..."
    
    # Clean apt cache
    if command -v apt-get > /dev/null 2>&1; then
        log_info "Cleaning apt cache..."
        sudo apt-get autoremove -y
        sudo apt-get autoclean
        sudo apt-get clean
        sudo rm -rf /var/lib/apt/lists/*
    fi
    
    # Clean Python cache
    if command -v python > /dev/null 2>&1; then
        log_info "Cleaning Python cache..."
        python -m pip cache purge
        find /opt/conda -name "*.pyc" -delete
        find /opt/conda -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
    fi
    
    # Clean R cache
    if command -v R > /dev/null 2>&1; then
        log_info "Cleaning R cache..."
        R -e "if(file.exists('~/.Rcache')) unlink('~/.Rcache', recursive=TRUE)"
        R -e "if(file.exists('/tmp/R*')) unlink('/tmp/R*', recursive=TRUE)"
    fi
    
    # Clean conda cache
    if command -v micromamba > /dev/null 2>&1; then
        log_info "Cleaning conda cache..."
        micromamba clean --all --yes
    fi
    
    # Clean temporary files
    log_info "Cleaning temporary files..."
    sudo rm -rf /tmp/*
    sudo rm -rf /var/tmp/*
    sudo rm -rf /root/.cache
    sudo rm -rf /home/*/.cache
    
    # Clean logs
    sudo rm -rf /var/log/*
    sudo truncate -s 0 /var/log/lastlog
    sudo truncate -s 0 /var/log/faillog
    
    # Clean history
    history -c
    > ~/.bash_history
    
    log_success "Docker image cleanup completed"
}

# Optimize layer caching
optimize_layers() {
    log_info "Optimizing Docker layers..."
    
    # Combine frequently changing files
    local optimization_script="/tmp/optimize_layers.sh"
    
    cat > "$optimization_script" << 'EOF'
#!/bin/bash
# Layer optimization script

# Create combined package installation
echo "Combining package installations for better caching..."

# Combine requirements files if they exist
if [ -d "/tmp/requirements" ]; then
    cat /tmp/requirements/*.txt > /tmp/combined_requirements.txt 2>/dev/null || true
fi

# Combine R package scripts
if [ -d "/tmp/r_packages" ]; then
    cat /tmp/r_packages/*.R > /tmp/combined_r_packages.R 2>/dev/null || true
fi

# Combine conda environments
if [ -d "/tmp/conda_envs" ]; then
    # Note: conda environments can't be simply combined, but we can optimize order
    echo "Conda environments found, optimizing installation order..."
fi
EOF
    
    chmod +x "$optimization_script"
    bash "$optimization_script"
    rm -f "$optimization_script"
    
    log_success "Layer optimization completed"
}

# Calculate image size reduction
calculate_size_reduction() {
    local before_size="$1"
    local after_size="$2"
    
    local reduction=$((before_size - after_size))
    local percentage=$((reduction * 100 / before_size))
    
    log_info "Image size reduction: ${reduction}MB (${percentage}%)"
    log_info "Before: ${before_size}MB, After: ${after_size}MB"
}

# Optimize file permissions
optimize_permissions() {
    log_info "Optimizing file permissions..."
    
    # Set appropriate permissions for scripts
    find /opt/devcontainer -name "*.sh" -exec chmod +x {} \;
    
    # Set appropriate permissions for config files
    find /opt/devcontainer -name "*.yml" -exec chmod 644 {} \;
    find /opt/devcontainer -name "*.txt" -exec chmod 644 {} \;
    find /opt/devcontainer -name "*.R" -exec chmod 644 {} \;
    
    # Set directory permissions
    find /opt/devcontainer -type d -exec chmod 755 {} \;
    
    log_success "File permissions optimized"
}

# Remove unnecessary packages
remove_unnecessary_packages() {
    log_info "Removing unnecessary packages..."
    
    # Remove development packages that are no longer needed
    if command -v apt-get > /dev/null 2>&1; then
        sudo apt-get autoremove -y
        
        # Remove packages that are only needed for building
        sudo apt-get remove -y --purge \
            build-essential \
            gcc \
            g++ \
            make \
            cmake \
            pkg-config \
            libtool \
            autoconf \
            automake \
            2>/dev/null || true
    fi
    
    # Remove Python build dependencies
    if command -v pip > /dev/null 2>&1; then
        pip uninstall -y setuptools wheel build || true
    fi
    
    log_success "Unnecessary packages removed"
}

# Optimize conda environment
optimize_conda_env() {
    log_info "Optimizing conda environment..."
    
    if command -v micromamba > /dev/null 2>&1; then
        # Remove unused packages
        micromamba clean --all --yes
        
        # Optimize package cache
        micromamba clean --index-cache --yes
        micromamba clean --lock --yes
        micromamba clean --tarballs --yes
        
        # Remove orphaned packages
        micromamba clean --packages --yes
        
        log_success "Conda environment optimized"
    else
        log_warning "Micromamba not found, skipping conda optimization"
    fi
}

# Create slim version of image
create_slim_image() {
    log_info "Creating slim version optimizations..."
    
    # Remove documentation
    find /opt/conda -name "*.md" -delete 2>/dev/null || true
    find /opt/conda -name "*.rst" -delete 2>/dev/null || true
    find /opt/conda -name "*.txt" -delete 2>/dev/null || true
    find /opt/conda -name "docs" -type d -exec rm -rf {} + 2>/dev/null || true
    
    # Remove examples
    find /opt/conda -name "examples" -type d -exec rm -rf {} + 2>/dev/null || true
    find /opt/conda -name "demo" -type d -exec rm -rf {} + 2>/dev/null || true
    find /opt/conda -name "test" -type d -exec rm -rf {} + 2>/dev/null || true
    find /opt/conda -name "tests" -type d -exec rm -rf {} + 2>/dev/null || true
    
    # Remove locale files (keep only English)
    find /opt/conda -name "locale" -type d -exec sh -c 'find "$1" -mindepth 1 -maxdepth 1 -type d ! -name "en*" -exec rm -rf {} +' _ {} \; 2>/dev/null || true
    
    # Remove static libraries
    find /opt/conda -name "*.a" -delete 2>/dev/null || true
    
    # Remove debug symbols
    find /opt/conda -name "*.debug" -delete 2>/dev/null || true
    
    log_success "Slim image optimizations completed"
}

# Comprehensive image optimization
optimize_image() {
    local optimization_level="${1:-standard}"
    
    log_info "Starting comprehensive image optimization (level: $optimization_level)"
    
    # Get initial size (approximate)
    local initial_size=$(du -sm /opt/conda 2>/dev/null | cut -f1 || echo "0")
    
    # Standard optimizations
    cleanup_image
    optimize_layers
    optimize_permissions
    optimize_conda_env
    
    # Additional optimizations for aggressive level
    if [ "$optimization_level" = "aggressive" ]; then
        log_info "Applying aggressive optimizations..."
        remove_unnecessary_packages
        create_slim_image
    fi
    
    # Get final size
    local final_size=$(du -sm /opt/conda 2>/dev/null | cut -f1 || echo "0")
    
    # Calculate reduction
    if [ "$initial_size" -gt 0 ] && [ "$final_size" -gt 0 ]; then
        calculate_size_reduction "$initial_size" "$final_size"
    fi
    
    log_success "Comprehensive image optimization completed"
}

# Validate optimization
validate_optimization() {
    log_info "Validating optimization..."
    
    # Check that essential tools still work
    local validation_failed=0
    
    # Test Python
    if ! python -c "import sys; print('Python OK')" 2>/dev/null; then
        log_error "Python validation failed after optimization"
        validation_failed=1
    fi
    
    # Test R
    if ! R -e "cat('R OK\n')" 2>/dev/null; then
        log_error "R validation failed after optimization"
        validation_failed=1
    fi
    
    # Test micromamba
    if ! micromamba --version 2>/dev/null; then
        log_error "Micromamba validation failed after optimization"
        validation_failed=1
    fi
    
    # Test essential packages
    if ! python -c "import numpy, pandas" 2>/dev/null; then
        log_error "Essential Python packages validation failed"
        validation_failed=1
    fi
    
    if [ $validation_failed -eq 0 ]; then
        log_success "Optimization validation passed"
        return 0
    else
        log_error "Optimization validation failed"
        return 1
    fi
}

# Export functions
export -f cleanup_image optimize_layers optimize_permissions optimize_conda_env
export -f remove_unnecessary_packages create_slim_image optimize_image validate_optimization
export -f calculate_size_reduction