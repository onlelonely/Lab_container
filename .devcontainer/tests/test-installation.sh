#!/bin/bash

# Installation test suite
set -euo pipefail

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UTILS_DIR="$SCRIPT_DIR/../scripts/utils"

# Source utility functions
source "$UTILS_DIR/logging.sh"
source "$UTILS_DIR/validation.sh"

# Script name for logging
readonly SCRIPT_NAME="test-installation.sh"

# Test configuration
readonly TEST_RESULTS_DIR="/workspace/logs/test-results"
readonly TEST_RESULTS_FILE="$TEST_RESULTS_DIR/installation-test-$(date +%Y%m%d-%H%M%S).log"

# Create test results directory
mkdir -p "$TEST_RESULTS_DIR"

# Test functions
test_python_packages() {
    log_info "Testing Python packages..."
    # Essential packages that should be available in core target (from requirements-core.txt)
    local essential_packages=()  # Actually core target has no data science packages
    # All data science packages are optional in core target
    local optional_packages=("numpy" "pandas" "matplotlib" "scipy" "scikit-learn" "jupyter" "biopython")
    local failed_core=()
    local failed_optional=()
    
    # Test essential packages (none in core target)
    for package in "${essential_packages[@]}"; do
        if python -c "import $package" 2>/dev/null; then
            log_success "✓ $package installed correctly"
            echo "PASS: Python package $package" >> "$TEST_RESULTS_FILE"
        else
            log_error "✗ $package installation failed"
            echo "FAIL: Python package $package" >> "$TEST_RESULTS_FILE"
            failed_core+=("$package")
        fi
    done
    
    # Test optional packages
    for package in "${optional_packages[@]}"; do
        if python -c "import $package" 2>/dev/null; then
            log_success "✓ $package installed correctly"
            echo "PASS: Python package $package" >> "$TEST_RESULTS_FILE"
        else
            log_warning "⚠ $package not available (optional for core target)"
            echo "WARN: Python package $package not installed" >> "$TEST_RESULTS_FILE"
            failed_optional+=("$package")
        fi
    done
    
    if [ ${#failed_core[@]} -eq 0 ]; then
        log_success "Python packages test passed (core target has minimal packages)"
        if [ ${#failed_optional[@]} -gt 0 ]; then
            log_info "Optional packages missing: ${failed_optional[*]} (expected for core target)"
        fi
        return 0
    else
        log_error "Failed essential Python packages: ${failed_core[*]}"
        return 1
    fi
}

test_r_packages() {
    log_info "Testing R packages..."
    # Core R packages that should be available in core target
    local core_packages=("base" "utils" "stats" "data.table" "devtools" "rmarkdown" "knitr")
    # Optional packages that may not be in core target
    local optional_packages=("tidyverse" "ggplot2" "Biostrings" "DESeq2" "dplyr")
    local failed_core=()
    local failed_optional=()
    
    # Test core R packages
    for package in "${core_packages[@]}"; do
        if R -e "library($package)" 2>/dev/null >/dev/null; then
            log_success "✓ $package installed correctly"
            echo "PASS: R package $package" >> "$TEST_RESULTS_FILE"
        else
            log_error "✗ $package installation failed"
            echo "FAIL: R package $package" >> "$TEST_RESULTS_FILE"
            failed_core+=("$package")
        fi
    done
    
    # Test optional R packages
    for package in "${optional_packages[@]}"; do
        if R -e "library($package)" 2>/dev/null >/dev/null; then
            log_success "✓ $package installed correctly"
            echo "PASS: R package $package" >> "$TEST_RESULTS_FILE"
        else
            log_warning "⚠ $package not available (optional for core target)"
            echo "WARN: R package $package not installed" >> "$TEST_RESULTS_FILE"
            failed_optional+=("$package")
        fi
    done
    
    if [ ${#failed_core[@]} -eq 0 ]; then
        log_success "Core R packages test passed"
        if [ ${#failed_optional[@]} -gt 0 ]; then
            log_info "Optional R packages missing: ${failed_optional[*]} (expected for core target)"
        fi
        return 0
    else
        log_error "Failed core R packages: ${failed_core[*]}"
        return 1
    fi
}

test_conda_packages() {
    log_info "Testing conda packages..."
    # Essential packages that must be available in core target (from environment-core.yml)
    local essential_packages=("python" "r-base" "jupyter" "pip" "git" "curl" "wget")
    # Optional packages that may not be in core target
    local optional_packages=()
    local failed_essential=()
    local failed_optional=()
    
    # Test essential packages
    for package in "${essential_packages[@]}"; do
        if conda list | grep -q "^$package"; then
            log_success "✓ $package installed correctly"
            echo "PASS: Conda package $package" >> "$TEST_RESULTS_FILE"
        else
            log_error "✗ $package installation failed"
            echo "FAIL: Conda package $package" >> "$TEST_RESULTS_FILE"
            failed_essential+=("$package")
        fi
    done
    
    # Test optional packages
    for package in "${optional_packages[@]}"; do
        if conda list | grep -q "^$package"; then
            log_success "✓ $package installed correctly"
            echo "PASS: Conda package $package" >> "$TEST_RESULTS_FILE"
        else
            log_warning "⚠ $package not available (optional for core target)"
            echo "WARN: Conda package $package not installed" >> "$TEST_RESULTS_FILE"
            failed_optional+=("$package")
        fi
    done
    
    if [ ${#failed_essential[@]} -eq 0 ]; then
        log_success "Essential conda packages test passed"
        if [ ${#optional_packages[@]} -gt 0 ] && [ ${#failed_optional[@]} -gt 0 ]; then
            log_info "Optional conda packages missing: ${failed_optional[*]} (expected for core target)"
        fi
        return 0
    else
        log_error "Failed essential conda packages: ${failed_essential[*]}"
        return 1
    fi
}

test_jupyter_functionality() {
    log_info "Testing Jupyter functionality..."
    
    # Check if Jupyter is available first
    if ! command -v jupyter > /dev/null 2>&1; then
        log_warning "⚠ Jupyter not installed (optional for core target)"
        echo "WARN: Jupyter not available" >> "$TEST_RESULTS_FILE"
        return 0
    fi
    
    # Test Jupyter kernels
    local kernels=$(jupyter kernelspec list 2>/dev/null | grep -E "(python|ir)" | wc -l)
    if [ "$kernels" -ge 1 ]; then
        log_success "✓ Jupyter kernels available: $kernels"
        echo "PASS: Jupyter kernels ($kernels available)" >> "$TEST_RESULTS_FILE"
    else
        log_warning "⚠ Limited Jupyter kernels: $kernels"
        echo "WARN: Jupyter kernels ($kernels available)" >> "$TEST_RESULTS_FILE"
    fi
    
    # Test Jupyter version (lighter test than full startup)
    if jupyter --version > /dev/null 2>&1; then
        log_success "✓ Jupyter is functional"
        echo "PASS: Jupyter functionality" >> "$TEST_RESULTS_FILE"
    else
        log_warning "⚠ Jupyter version check failed"
        echo "WARN: Jupyter functionality" >> "$TEST_RESULTS_FILE"
    fi
    
    return 0
}

test_development_tools() {
    log_info "Testing development tools..."
    local tools=("git" "curl" "wget" "python" "R" "conda")
    local failed_tools=()
    
    for tool in "${tools[@]}"; do
        if command -v "$tool" > /dev/null 2>&1; then
            log_success "✓ $tool is available"
            echo "PASS: Development tool $tool" >> "$TEST_RESULTS_FILE"
        else
            log_error "✗ $tool is not available"
            echo "FAIL: Development tool $tool" >> "$TEST_RESULTS_FILE"
            failed_tools+=("$tool")
        fi
    done
    
    if [ ${#failed_tools[@]} -eq 0 ]; then
        log_success "All development tools test passed"
        return 0
    else
        log_error "Failed development tools: ${failed_tools[*]}"
        return 1
    fi
}

test_file_system() {
    log_info "Testing file system..."
    local directories=("/workspace" "/workspace/logs" "/opt/conda")
    local failed_directories=()
    
    for dir in "${directories[@]}"; do
        if [ -d "$dir" ] && [ -w "$dir" ]; then
            log_success "✓ Directory accessible: $dir"
            echo "PASS: Directory $dir" >> "$TEST_RESULTS_FILE"
        else
            log_error "✗ Directory not accessible: $dir"
            echo "FAIL: Directory $dir" >> "$TEST_RESULTS_FILE"
            failed_directories+=("$dir")
        fi
    done
    
    if [ ${#failed_directories[@]} -eq 0 ]; then
        log_success "All file system tests passed"
        return 0
    else
        log_error "Failed directory tests: ${failed_directories[*]}"
        return 1
    fi
}

# Performance test
test_performance() {
    log_info "Testing performance..."
    
    # Test Python import speed
    local python_time=$(time (python -c "import numpy, pandas, matplotlib" 2>/dev/null) 2>&1 | grep real | cut -d' ' -f2)
    log_info "Python import time: $python_time"
    echo "INFO: Python import time: $python_time" >> "$TEST_RESULTS_FILE"
    
    # Test R library load speed
    local r_time=$(time (R -e "library(tidyverse); library(ggplot2)" 2>/dev/null) 2>&1 | grep real | cut -d' ' -f2)
    log_info "R library load time: $r_time"
    echo "INFO: R library load time: $r_time" >> "$TEST_RESULTS_FILE"
    
    return 0
}

# Generate test report
generate_test_report() {
    local total_tests=$(grep -c "PASS\|FAIL" "$TEST_RESULTS_FILE")
    local passed_tests=$(grep -c "PASS" "$TEST_RESULTS_FILE")
    local failed_tests=$(grep -c "FAIL" "$TEST_RESULTS_FILE")
    
    log_info "Test Results Summary:"
    log_info "Total tests: $total_tests"
    log_info "Passed: $passed_tests"
    log_info "Failed: $failed_tests"
    
    echo "=== TEST SUMMARY ===" >> "$TEST_RESULTS_FILE"
    echo "Total tests: $total_tests" >> "$TEST_RESULTS_FILE"
    echo "Passed: $passed_tests" >> "$TEST_RESULTS_FILE"
    echo "Failed: $failed_tests" >> "$TEST_RESULTS_FILE"
    echo "Test completed at: $(date)" >> "$TEST_RESULTS_FILE"
    
    if [ "$failed_tests" -eq 0 ]; then
        log_success "All tests passed!"
        return 0
    else
        log_error "Some tests failed. Check $TEST_RESULTS_FILE for details."
        return 1
    fi
}

# Main test function
main() {
    log_script_start "$SCRIPT_NAME"
    
    echo "=== INSTALLATION TEST REPORT ===" > "$TEST_RESULTS_FILE"
    echo "Test started at: $(date)" >> "$TEST_RESULTS_FILE"
    echo "" >> "$TEST_RESULTS_FILE"
    
    local test_failed=0
    
    # Run all tests
    test_file_system || test_failed=1
    test_development_tools || test_failed=1
    test_conda_packages || test_failed=1
    test_python_packages || test_failed=1
    test_r_packages || test_failed=1
    test_jupyter_functionality || test_failed=1
    test_performance
    
    # Generate report
    generate_test_report
    
    log_script_end "$SCRIPT_NAME" $test_failed
    
    return $test_failed
}

# Error handling
trap 'log_script_end "$SCRIPT_NAME" $?' EXIT

# Run main function
main "$@"