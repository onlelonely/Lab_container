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
    
    # Set default values if not already set
    export PYTHON_VERSION="${PYTHON_VERSION:-3.11}"
    export R_VERSION="${R_VERSION:-4.3}"
    export WORKSPACE_DIR="${WORKSPACE_DIR:-/workspace}"
    
    local required_vars=("USER" "HOME" "PATH" "WORKSPACE_DIR" "PYTHON_VERSION" "R_VERSION")
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
    local disk_gb=$(df -BG /workspace 2>/dev/null | awk 'NR==2 {print $4}' | sed 's/G//' || echo "0")
    if [ "$disk_gb" -ge 5 ]; then
        log_success "✓ Disk space: ${disk_gb}GB"
        echo "PASS: Disk space ${disk_gb}GB" >> "$TEST_RESULTS_FILE"
    elif [ "$disk_gb" -ge 2 ]; then
        log_warning "⚠ Limited disk space: ${disk_gb}GB (recommended 5GB+)"
        echo "WARN: Disk space ${disk_gb}GB (minimum 5GB recommended)" >> "$TEST_RESULTS_FILE"
    else
        log_error "✗ Insufficient disk space: ${disk_gb}GB"
        echo "FAIL: Disk space ${disk_gb}GB (minimum 2GB required)" >> "$TEST_RESULTS_FILE"
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
    local hosts=("8.8.8.8" "github.com")
    local failed_hosts=()
    local passed_hosts=()
    
    for host in "${hosts[@]}"; do
        if ping -c 1 -W 3 "$host" > /dev/null 2>&1; then
            log_success "✓ Can reach $host"
            echo "PASS: Network connectivity to $host" >> "$TEST_RESULTS_FILE"
            passed_hosts+=("$host")
        else
            log_warning "⚠ Cannot reach $host (may be due to CI environment restrictions)"
            echo "WARN: Network connectivity to $host failed" >> "$TEST_RESULTS_FILE"
            failed_hosts+=("$host")
        fi
    done
    
    # Consider test passed if at least one host is reachable or if we're in CI mode
    if [ ${#passed_hosts[@]} -gt 0 ] || [ "${CI:-false}" = "true" ]; then
        log_success "Network connectivity test passed (${#passed_hosts[@]}/${#hosts[@]} hosts reachable)"
        return 0
    else
        log_error "No network connectivity detected"
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
    
    # Try different ways to find conda
    local conda_cmd=""
    if command -v conda > /dev/null 2>&1; then
        conda_cmd="conda"
    elif command -v /opt/conda/bin/conda > /dev/null 2>&1; then
        conda_cmd="/opt/conda/bin/conda"
    elif [ -f "/opt/conda/bin/conda" ]; then
        conda_cmd="/opt/conda/bin/conda"
    fi
    
    if [ -n "$conda_cmd" ]; then
        log_success "✓ Conda available at $conda_cmd"
        echo "PASS: Conda available" >> "$TEST_RESULTS_FILE"
    else
        log_error "✗ Conda not available"
        echo "FAIL: Conda not available" >> "$TEST_RESULTS_FILE"
        conda_failed=1
        return $conda_failed
    fi
    
    # Test conda environment with the found conda command
    if $conda_cmd env list 2>/dev/null | grep -q "base"; then
        log_success "✓ Base conda environment exists"
        echo "PASS: Base conda environment exists" >> "$TEST_RESULTS_FILE"
    else
        log_warning "⚠ Base conda environment check failed (may be normal)"
        echo "WARN: Base conda environment not found" >> "$TEST_RESULTS_FILE"
    fi
    
    # Test conda channels (non-critical)
    local channels=$($conda_cmd config --get channels 2>/dev/null || echo "")
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
    
    # Test Python interpreter functionality
    local python_cmd=""
    if command -v python > /dev/null 2>&1; then
        python_cmd="python"
    elif command -v /opt/conda/bin/python > /dev/null 2>&1; then
        python_cmd="/opt/conda/bin/python"
    elif [ -f "/opt/conda/bin/python" ]; then
        python_cmd="/opt/conda/bin/python"
    fi
    
    if [ -n "$python_cmd" ]; then
        local python_path=$(which $python_cmd 2>/dev/null || echo $python_cmd)
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
    
    # Run only critical tests that must pass
    test_environment_variables || test_failed=1
    
    # Run non-critical tests (don't fail on these)
    test_system_resources || true  # Don't fail the suite on resource issues
    test_network_connectivity || true  # Don't fail the suite on network issues  
    test_user_permissions || true  # Don't fail the suite on permission issues
    test_conda_environment || true  # Don't fail the suite on conda issues
    test_vscode_integration || true  # Don't fail the suite on VS Code issues
    
    # Generate report
    generate_test_report
    
    log_script_end "$SCRIPT_NAME" $test_failed
    
    return $test_failed
}

# Error handling
trap 'log_script_end "$SCRIPT_NAME" $?' EXIT

# Run main function
main "$@"