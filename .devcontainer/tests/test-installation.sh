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
    local packages=("numpy" "pandas" "matplotlib" "scipy" "scikit-learn" "jupyter" "biopython")
    local failed_packages=()
    
    for package in "${packages[@]}"; do
        if python -c "import $package" 2>/dev/null; then
            log_success "✓ $package installed correctly"
            echo "PASS: Python package $package" >> "$TEST_RESULTS_FILE"
        else
            log_error "✗ $package installation failed"
            echo "FAIL: Python package $package" >> "$TEST_RESULTS_FILE"
            failed_packages+=("$package")
        fi
    done
    
    if [ ${#failed_packages[@]} -eq 0 ]; then
        log_success "All Python packages test passed"
        return 0
    else
        log_error "Failed Python packages: ${failed_packages[*]}"
        return 1
    fi
}

test_r_packages() {
    log_info "Testing R packages..."
    local packages=("tidyverse" "ggplot2" "data.table" "Biostrings" "DESeq2" "dplyr")
    local failed_packages=()
    
    for package in "${packages[@]}"; do
        if R -e "library($package)" 2>/dev/null; then
            log_success "✓ $package installed correctly"
            echo "PASS: R package $package" >> "$TEST_RESULTS_FILE"
        else
            log_error "✗ $package installation failed"
            echo "FAIL: R package $package" >> "$TEST_RESULTS_FILE"
            failed_packages+=("$package")
        fi
    done
    
    if [ ${#failed_packages[@]} -eq 0 ]; then
        log_success "All R packages test passed"
        return 0
    else
        log_error "Failed R packages: ${failed_packages[*]}"
        return 1
    fi
}

test_conda_packages() {
    log_info "Testing conda packages..."
    local packages=("python" "r-base" "jupyter" "git" "curl" "wget")
    local failed_packages=()
    
    for package in "${packages[@]}"; do
        if conda list | grep -q "^$package"; then
            log_success "✓ $package installed correctly"
            echo "PASS: Conda package $package" >> "$TEST_RESULTS_FILE"
        else
            log_error "✗ $package installation failed"
            echo "FAIL: Conda package $package" >> "$TEST_RESULTS_FILE"
            failed_packages+=("$package")
        fi
    done
    
    if [ ${#failed_packages[@]} -eq 0 ]; then
        log_success "All conda packages test passed"
        return 0
    else
        log_error "Failed conda packages: ${failed_packages[*]}"
        return 1
    fi
}

test_jupyter_functionality() {
    log_info "Testing Jupyter functionality..."
    
    # Test Jupyter kernels
    local kernels=$(jupyter kernelspec list 2>/dev/null | grep -E "(python|ir)" | wc -l)
    if [ "$kernels" -ge 2 ]; then
        log_success "✓ Jupyter kernels available: $kernels"
        echo "PASS: Jupyter kernels ($kernels available)" >> "$TEST_RESULTS_FILE"
    else
        log_error "✗ Insufficient Jupyter kernels: $kernels"
        echo "FAIL: Jupyter kernels ($kernels available)" >> "$TEST_RESULTS_FILE"
        return 1
    fi
    
    # Test Jupyter startup
    if timeout 30 jupyter lab --version > /dev/null 2>&1; then
        log_success "✓ Jupyter Lab is functional"
        echo "PASS: Jupyter Lab functionality" >> "$TEST_RESULTS_FILE"
    else
        log_error "✗ Jupyter Lab is not functional"
        echo "FAIL: Jupyter Lab functionality" >> "$TEST_RESULTS_FILE"
        return 1
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