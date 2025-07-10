#!/bin/bash

# Environment test suite
set -euo pipefail

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UTILS_DIR="$SCRIPT_DIR/../scripts/utils"

# Source utility functions
source "$UTILS_DIR/logging.sh"
source "$UTILS_DIR/validation.sh"

# Script name for logging
readonly SCRIPT_NAME="test-environment.sh"

# Test configuration
readonly TEST_RESULTS_DIR="/workspace/logs/test-results"
readonly TEST_RESULTS_FILE="$TEST_RESULTS_DIR/environment-test-$(date +%Y%m%d-%H%M%S).log"

# Create test results directory
mkdir -p "$TEST_RESULTS_DIR"

# Test environment variables
test_environment_variables() {
    log_info "Testing environment variables..."
    local required_vars=("USER" "HOME" "PATH" "WORKSPACE_DIR")
    local failed_vars=()
    
    for var in "${required_vars[@]}"; do
        if [ -n "${!var:-}" ]; then
            log_success "✓ $var = ${!var}"
            echo "PASS: Environment variable $var = ${!var}" >> "$TEST_RESULTS_FILE"
        else
            log_error "✗ $var is not set"
            echo "FAIL: Environment variable $var is not set" >> "$TEST_RESULTS_FILE"
            failed_vars+=("$var")
        fi
    done
    
    if [ ${#failed_vars[@]} -eq 0 ]; then
        log_success "All environment variables test passed"
        return 0
    else
        log_error "Failed environment variables: ${failed_vars[*]}"
        return 1
    fi
}

# Test system resources
test_system_resources() {
    log_info "Testing system resources..."
    local resource_failed=0
    
    # Test memory
    local memory_gb=$(free -g | awk '/^Mem:/{print $2}')
    if [ "$memory_gb" -ge 2 ]; then
        log_success "✓ Memory: ${memory_gb}GB"
        echo "PASS: Memory ${memory_gb}GB" >> "$TEST_RESULTS_FILE"
    else
        log_error "✗ Insufficient memory: ${memory_gb}GB"
        echo "FAIL: Memory ${memory_gb}GB (minimum 2GB required)" >> "$TEST_RESULTS_FILE"
        resource_failed=1
    fi
    
    # Test disk space
    local disk_gb=$(df -BG /workspace | awk 'NR==2 {print $4}' | sed 's/G//')
    if [ "$disk_gb" -ge 5 ]; then
        log_success "✓ Disk space: ${disk_gb}GB"
        echo "PASS: Disk space ${disk_gb}GB" >> "$TEST_RESULTS_FILE"
    else
        log_error "✗ Insufficient disk space: ${disk_gb}GB"
        echo "FAIL: Disk space ${disk_gb}GB (minimum 5GB required)" >> "$TEST_RESULTS_FILE"
        resource_failed=1
    fi
    
    # Test CPU cores
    local cpu_cores=$(nproc)
    if [ "$cpu_cores" -ge 2 ]; then
        log_success "✓ CPU cores: $cpu_cores"
        echo "PASS: CPU cores $cpu_cores" >> "$TEST_RESULTS_FILE"
    else
        log_warning "⚠ Limited CPU cores: $cpu_cores"
        echo "WARN: CPU cores $cpu_cores (recommended 2+ cores)" >> "$TEST_RESULTS_FILE"
    fi
    
    return $resource_failed
}

# Test network connectivity
test_network_connectivity() {
    log_info "Testing network connectivity..."
    local hosts=("8.8.8.8" "github.com" "conda-forge.org" "pypi.org")
    local failed_hosts=()
    
    for host in "${hosts[@]}"; do
        if ping -c 1 -W 5 "$host" > /dev/null 2>&1; then
            log_success "✓ Can reach $host"
            echo "PASS: Network connectivity to $host" >> "$TEST_RESULTS_FILE"
        else
            log_error "✗ Cannot reach $host"
            echo "FAIL: Network connectivity to $host" >> "$TEST_RESULTS_FILE"
            failed_hosts+=("$host")
        fi
    done
    
    if [ ${#failed_hosts[@]} -eq 0 ]; then
        log_success "All network connectivity tests passed"
        return 0
    else
        log_error "Failed network connectivity to: ${failed_hosts[*]}"
        return 1
    fi
}

# Test user permissions
test_user_permissions() {
    log_info "Testing user permissions..."
    local permission_failed=0
    
    # Test workspace write permissions
    local test_file="/workspace/.test_write_permission"
    if touch "$test_file" 2>/dev/null; then
        log_success "✓ Workspace write permission"
        echo "PASS: Workspace write permission" >> "$TEST_RESULTS_FILE"
        rm -f "$test_file"
    else
        log_error "✗ No workspace write permission"
        echo "FAIL: Workspace write permission" >> "$TEST_RESULTS_FILE"
        permission_failed=1
    fi
    
    # Test user is not root
    if [ "$USER" != "root" ]; then
        log_success "✓ Non-root user: $USER"
        echo "PASS: Non-root user $USER" >> "$TEST_RESULTS_FILE"
    else
        log_error "✗ Running as root user"
        echo "FAIL: Running as root user" >> "$TEST_RESULTS_FILE"
        permission_failed=1
    fi
    
    # Test sudo access
    if sudo -n true 2>/dev/null; then
        log_success "✓ Sudo access available"
        echo "PASS: Sudo access available" >> "$TEST_RESULTS_FILE"
    else
        log_warning "⚠ No sudo access"
        echo "WARN: No sudo access" >> "$TEST_RESULTS_FILE"
    fi
    
    return $permission_failed
}

# Test conda environment
test_conda_environment() {
    log_info "Testing conda environment..."
    local conda_failed=0
    
    # Test micromamba command
    if command -v micromamba > /dev/null 2>&1; then
        log_success "✓ Micromamba available"
        echo "PASS: Micromamba available" >> "$TEST_RESULTS_FILE"
    else
        log_error "✗ Micromamba not available"
        echo "FAIL: Micromamba not available" >> "$TEST_RESULTS_FILE"
        conda_failed=1
    fi
    
    # Test conda environment
    if micromamba env list | grep -q "base"; then
        log_success "✓ Base conda environment exists"
        echo "PASS: Base conda environment exists" >> "$TEST_RESULTS_FILE"
    else
        log_error "✗ Base conda environment missing"
        echo "FAIL: Base conda environment missing" >> "$TEST_RESULTS_FILE"
        conda_failed=1
    fi
    
    # Test conda channels
    local channels=$(micromamba config --get channels 2>/dev/null || echo "")
    if echo "$channels" | grep -q "conda-forge"; then
        log_success "✓ Conda-forge channel configured"
        echo "PASS: Conda-forge channel configured" >> "$TEST_RESULTS_FILE"
    else
        log_warning "⚠ Conda-forge channel not configured"
        echo "WARN: Conda-forge channel not configured" >> "$TEST_RESULTS_FILE"
    fi
    
    return $conda_failed
}

# Test VS Code integration
test_vscode_integration() {
    log_info "Testing VS Code integration..."
    local vscode_failed=0
    
    # Test VS Code extensions directory
    local extensions_dir="/home/$USER/.vscode-server/extensions"
    if [ -d "$extensions_dir" ]; then
        log_success "✓ VS Code extensions directory exists"
        echo "PASS: VS Code extensions directory exists" >> "$TEST_RESULTS_FILE"
        
        # Count extensions
        local extension_count=$(find "$extensions_dir" -maxdepth 1 -type d | wc -l)
        log_info "VS Code extensions installed: $extension_count"
        echo "INFO: VS Code extensions installed: $extension_count" >> "$TEST_RESULTS_FILE"
    else
        log_warning "⚠ VS Code extensions directory not found"
        echo "WARN: VS Code extensions directory not found" >> "$TEST_RESULTS_FILE"
    fi
    
    # Test Python extension functionality
    if command -v python > /dev/null 2>&1; then
        local python_path=$(which python)
        log_success "✓ Python interpreter available at $python_path"
        echo "PASS: Python interpreter available at $python_path" >> "$TEST_RESULTS_FILE"
    else
        log_error "✗ Python interpreter not found"
        echo "FAIL: Python interpreter not found" >> "$TEST_RESULTS_FILE"
        vscode_failed=1
    fi
    
    return $vscode_failed
}

# Generate test report
generate_test_report() {
    local total_tests=$(grep -c "PASS\|FAIL" "$TEST_RESULTS_FILE")
    local passed_tests=$(grep -c "PASS" "$TEST_RESULTS_FILE")
    local failed_tests=$(grep -c "FAIL" "$TEST_RESULTS_FILE")
    local warnings=$(grep -c "WARN" "$TEST_RESULTS_FILE")
    
    log_info "Environment Test Results Summary:"
    log_info "Total tests: $total_tests"
    log_info "Passed: $passed_tests"
    log_info "Failed: $failed_tests"
    log_info "Warnings: $warnings"
    
    echo "=== ENVIRONMENT TEST SUMMARY ===" >> "$TEST_RESULTS_FILE"
    echo "Total tests: $total_tests" >> "$TEST_RESULTS_FILE"
    echo "Passed: $passed_tests" >> "$TEST_RESULTS_FILE"
    echo "Failed: $failed_tests" >> "$TEST_RESULTS_FILE"
    echo "Warnings: $warnings" >> "$TEST_RESULTS_FILE"
    echo "Test completed at: $(date)" >> "$TEST_RESULTS_FILE"
    
    if [ "$failed_tests" -eq 0 ]; then
        log_success "All environment tests passed!"
        return 0
    else
        log_error "Some environment tests failed. Check $TEST_RESULTS_FILE for details."
        return 1
    fi
}

# Main test function
main() {
    log_script_start "$SCRIPT_NAME"
    
    echo "=== ENVIRONMENT TEST REPORT ===" > "$TEST_RESULTS_FILE"
    echo "Test started at: $(date)" >> "$TEST_RESULTS_FILE"
    echo "" >> "$TEST_RESULTS_FILE"
    
    local test_failed=0
    
    # Run all tests
    test_environment_variables || test_failed=1
    test_system_resources || test_failed=1
    test_network_connectivity || test_failed=1
    test_user_permissions || test_failed=1
    test_conda_environment || test_failed=1
    test_vscode_integration || test_failed=1
    
    # Generate report
    generate_test_report
    
    log_script_end "$SCRIPT_NAME" $test_failed
    
    return $test_failed
}

# Error handling
trap 'log_script_end "$SCRIPT_NAME" $?' EXIT

# Run main function
main "$@"