#!/bin/bash

# Docker image analysis utilities
set -euo pipefail

# Source logging functions
source "$(dirname "${BASH_SOURCE[0]}")/logging.sh"

# Analyze image layers
analyze_layers() {
    local image_name="${1:-devcontainer}"
    
    log_info "Analyzing Docker image layers for $image_name"
    
    # Get image history
    if command -v docker > /dev/null 2>&1; then
        docker history "$image_name" --format "table {{.CreatedBy}}\t{{.Size}}" | head -20
    else
        log_warning "Docker not available for layer analysis"
    fi
}

# Analyze directory sizes
analyze_directory_sizes() {
    log_info "Analyzing directory sizes..."
    
    # Analyze major directories
    local dirs=(
        "/opt/conda"
        "/usr/local"
        "/var/lib"
        "/home"
        "/tmp"
    )
    
    for dir in "${dirs[@]}"; do
        if [ -d "$dir" ]; then
            local size=$(du -sh "$dir" 2>/dev/null | cut -f1)
            log_info "Directory $dir: $size"
        fi
    done
}

# Analyze package sizes
analyze_package_sizes() {
    log_info "Analyzing package sizes..."
    
    # Python packages
    if command -v pip > /dev/null 2>&1; then
        log_info "Top 10 largest Python packages:"
        pip list --format=columns | tail -n +3 | while read -r package version; do
            if [ -n "$package" ]; then
                local size=$(pip show "$package" 2>/dev/null | grep -i "location" | head -1 | cut -d' ' -f2)
                if [ -n "$size" ] && [ -d "$size" ]; then
                    local pkg_size=$(du -sh "$size" 2>/dev/null | cut -f1)
                    echo "$package: $pkg_size"
                fi
            fi
        done | sort -hr | head -10
    fi
    
    # Conda packages
    if command -v micromamba > /dev/null 2>&1; then
        log_info "Top 10 largest conda packages:"
        micromamba list --json 2>/dev/null | jq -r '.[] | "\(.name) \(.version)"' | while read -r package version; do
            local pkg_dir="/opt/conda/pkgs/${package}-${version}"
            if [ -d "$pkg_dir" ]; then
                local size=$(du -sh "$pkg_dir" 2>/dev/null | cut -f1)
                echo "$package: $size"
            fi
        done | sort -hr | head -10 2>/dev/null || log_warning "jq not available for conda analysis"
    fi
}

# Find large files
find_large_files() {
    local min_size="${1:-100M}"
    local search_path="${2:-/opt/conda}"
    
    log_info "Finding files larger than $min_size in $search_path"
    
    find "$search_path" -type f -size "+$min_size" -exec ls -lh {} \; 2>/dev/null | \
        awk '{print $5 " " $9}' | sort -hr | head -20
}

# Analyze cache sizes
analyze_cache_sizes() {
    log_info "Analyzing cache sizes..."
    
    local cache_dirs=(
        "/home/*/.cache"
        "/root/.cache"
        "/tmp"
        "/var/tmp"
        "/opt/conda/pkgs"
        "/opt/conda/conda-meta"
    )
    
    for cache_dir in "${cache_dirs[@]}"; do
        if ls -d $cache_dir 2>/dev/null; then
            local size=$(du -sh $cache_dir 2>/dev/null | cut -f1)
            log_info "Cache $cache_dir: $size"
        fi
    done
}

# Generate optimization report
generate_optimization_report() {
    local report_file="${1:-/workspace/logs/optimization-report.txt}"
    
    log_info "Generating optimization report: $report_file"
    
    mkdir -p "$(dirname "$report_file")"
    
    cat > "$report_file" << EOF
=== Docker Image Optimization Report ===
Generated: $(date)
System: $(uname -a)

=== Directory Sizes ===
EOF
    
    # Append directory analysis
    analyze_directory_sizes >> "$report_file" 2>&1
    
    cat >> "$report_file" << EOF

=== Large Files (>100MB) ===
EOF
    
    # Append large files
    find_large_files >> "$report_file" 2>&1
    
    cat >> "$report_file" << EOF

=== Cache Analysis ===
EOF
    
    # Append cache analysis
    analyze_cache_sizes >> "$report_file" 2>&1
    
    cat >> "$report_file" << EOF

=== Package Sizes ===
EOF
    
    # Append package analysis
    analyze_package_sizes >> "$report_file" 2>&1
    
    cat >> "$report_file" << EOF

=== Optimization Recommendations ===
EOF
    
    # Generate recommendations
    generate_recommendations >> "$report_file" 2>&1
    
    log_success "Optimization report generated: $report_file"
}

# Generate optimization recommendations
generate_recommendations() {
    log_info "Generating optimization recommendations..."
    
    # Check for common optimization opportunities
    local recommendations=()
    
    # Check for large cache directories
    if [ -d "/home/*/.cache" ] || [ -d "/root/.cache" ]; then
        recommendations+=("Clear user cache directories to save space")
    fi
    
    # Check for large package cache
    if [ -d "/opt/conda/pkgs" ]; then
        local pkg_size=$(du -sm "/opt/conda/pkgs" 2>/dev/null | cut -f1)
        if [ "$pkg_size" -gt 500 ]; then
            recommendations+=("Clean conda package cache (${pkg_size}MB)")
        fi
    fi
    
    # Check for development packages
    if command -v gcc > /dev/null 2>&1; then
        recommendations+=("Remove development packages (gcc, g++, make) if not needed")
    fi
    
    # Check for documentation
    if find /opt/conda -name "docs" -type d | head -1; then
        recommendations+=("Remove documentation directories to save space")
    fi
    
    # Check for test directories
    if find /opt/conda -name "test*" -type d | head -1; then
        recommendations+=("Remove test directories to save space")
    fi
    
    # Print recommendations
    if [ ${#recommendations[@]} -gt 0 ]; then
        for rec in "${recommendations[@]}"; do
            echo "- $rec"
        done
    else
        echo "No major optimization opportunities found"
    fi
}

# Compare image sizes
compare_image_sizes() {
    local before_image="$1"
    local after_image="$2"
    
    log_info "Comparing image sizes..."
    
    if command -v docker > /dev/null 2>&1; then
        local before_size=$(docker image inspect "$before_image" --format='{{.Size}}' 2>/dev/null || echo "0")
        local after_size=$(docker image inspect "$after_image" --format='{{.Size}}' 2>/dev/null || echo "0")
        
        if [ "$before_size" -gt 0 ] && [ "$after_size" -gt 0 ]; then
            local reduction=$((before_size - after_size))
            local percentage=$((reduction * 100 / before_size))
            
            log_info "Before: $before_image = $(numfmt --to=si $before_size)B"
            log_info "After: $after_image = $(numfmt --to=si $after_size)B"
            log_info "Reduction: $(numfmt --to=si $reduction)B ($percentage%)"
        else
            log_warning "Unable to get image sizes for comparison"
        fi
    else
        log_warning "Docker not available for size comparison"
    fi
}

# Main analysis function
analyze_image() {
    local analysis_type="${1:-full}"
    local report_file="${2:-/workspace/logs/image-analysis.txt}"
    
    log_info "Starting image analysis (type: $analysis_type)"
    
    case "$analysis_type" in
        "layers")
            analyze_layers
            ;;
        "sizes")
            analyze_directory_sizes
            ;;
        "packages")
            analyze_package_sizes
            ;;
        "cache")
            analyze_cache_sizes
            ;;
        "large-files")
            find_large_files
            ;;
        "full")
            generate_optimization_report "$report_file"
            ;;
        *)
            log_error "Unknown analysis type: $analysis_type"
            log_info "Available types: layers, sizes, packages, cache, large-files, full"
            return 1
            ;;
    esac
    
    log_success "Image analysis completed"
}

# Export functions
export -f analyze_layers analyze_directory_sizes analyze_package_sizes
export -f find_large_files analyze_cache_sizes generate_optimization_report
export -f generate_recommendations compare_image_sizes analyze_image