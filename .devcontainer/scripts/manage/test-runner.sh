#!/bin/bash

# Automated test runner for dev container
set -euo pipefail

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UTILS_DIR="$SCRIPT_DIR/../utils"
TESTS_DIR="$SCRIPT_DIR/../../tests"

# Source utility functions
source "$UTILS_DIR/logging.sh"
source "$UTILS_DIR/validation.sh"

# Script name for logging
readonly SCRIPT_NAME="test-runner.sh"

# Test configuration
readonly TEST_RESULTS_DIR="/workspace/logs/test-results"
readonly TEST_SUITE_FILE="$TEST_RESULTS_DIR/test-suite-$(date +%Y%m%d-%H%M%S).log"
readonly CI_MODE="${CI:-false}"

# Create test results directory
mkdir -p "$TEST_RESULTS_DIR"

# Set default environment variables if not already set
export PYTHON_VERSION="${PYTHON_VERSION:-3.11}"
export R_VERSION="${R_VERSION:-4.3}"
export WORKSPACE_DIR="${WORKSPACE_DIR:-/workspace}"

# Test suite configuration
declare -A TEST_SUITES=(
    ["quick"]="test-environment.sh"
    ["installation"]="test-installation.sh test-environment.sh"
    ["full"]="test-environment.sh test-installation.sh"
    ["performance"]="test-environment.sh test-installation.sh"
    ["ci"]="test-environment.sh test-installation.sh"
)

# Main test runner function
main() {
    log_script_start "$SCRIPT_NAME"
    
    # Parse command line arguments
    local test_suite="${1:-quick}"
    local verbose="${2:-false}"
    local fail_fast="${3:-false}"
    
    log_info "Running test suite: $test_suite"
    log_info "Verbose mode: $verbose"
    log_info "Fail fast mode: $fail_fast"
    
    # Initialize test suite log
    echo "=== AUTOMATED TEST SUITE REPORT ===" > "$TEST_SUITE_FILE"
    echo "Suite: $test_suite" >> "$TEST_SUITE_FILE"
    echo "Started: $(date)" >> "$TEST_SUITE_FILE"
    echo "CI Mode: $CI_MODE" >> "$TEST_SUITE_FILE"
    echo "" >> "$TEST_SUITE_FILE"
    
    # Validate test suite
    if [[ ! ${TEST_SUITES[$test_suite]+_} ]]; then
        log_error "Unknown test suite: $test_suite"
        log_info "Available suites: ${!TEST_SUITES[@]}"
        exit 1
    fi
    
    # Run pre-test validation
    if ! run_pre_test_validation; then
        log_error "Pre-test validation failed"
        exit 1
    fi
    
    # Run test suite
    local test_files=(${TEST_SUITES[$test_suite]})
    local total_tests=${#test_files[@]}
    local passed_tests=0
    local failed_tests=0
    
    log_info "Running $total_tests test files..."
    
    for test_file in "${test_files[@]}"; do
        local test_path="$TESTS_DIR/$test_file"
        
        if [ ! -f "$test_path" ]; then
            log_error "Test file not found: $test_path"
            failed_tests=$((failed_tests + 1))
            continue
        fi
        
        log_info "Running test: $test_file"
        
        # Run individual test
        if run_individual_test "$test_path" "$verbose"; then
            passed_tests=$((passed_tests + 1))
            log_success "Test passed: $test_file"
        else
            failed_tests=$((failed_tests + 1))
            log_error "Test failed: $test_file"
            
            if [ "$fail_fast" = "true" ]; then
                log_error "Fail fast mode enabled, stopping test suite"
                break
            fi
        fi
    done
    
    # Run post-test analysis
    run_post_test_analysis "$test_suite"
    
    # Generate final report
    generate_final_report "$total_tests" "$passed_tests" "$failed_tests"
    
    # Determine exit code
    if [ $failed_tests -eq 0 ]; then
        log_success "All tests passed!"
        log_script_end "$SCRIPT_NAME" 0
        exit 0
    else
        log_error "$failed_tests out of $total_tests tests failed"
        log_script_end "$SCRIPT_NAME" 1
        exit 1
    fi
}

# Run pre-test validation
run_pre_test_validation() {
    log_info "Running pre-test validation..."
    
    # Check environment
    if ! validate_environment; then
        log_error "Environment validation failed"
        return 1
    fi
    
    # Check system resources
    validate_system_resources
    
    # Check test files exist
    for suite in "${!TEST_SUITES[@]}"; do
        local test_files=(${TEST_SUITES[$suite]})
        for test_file in "${test_files[@]}"; do
            if [ ! -f "$TESTS_DIR/$test_file" ]; then
                log_error "Test file missing: $TESTS_DIR/$test_file"
                return 1
            fi
        done
    done
    
    log_success "Pre-test validation passed"
    return 0
}

# Run individual test
run_individual_test() {
    local test_path="$1"
    local verbose="$2"
    local test_name=$(basename "$test_path")
    
    # Create individual test log
    local test_log="$TEST_RESULTS_DIR/${test_name%.sh}-$(date +%Y%m%d-%H%M%S).log"
    
    # Set up test environment
    export TEST_MODE="automated"
    export TEST_VERBOSE="$verbose"
    export TEST_LOG_FILE="$test_log"
    
    # Run test with timeout
    local timeout_duration=600  # 10 minutes
    local test_start_time=$(date +%s)
    
    if [ "$verbose" = "true" ]; then
        log_info "Running test with verbose output: $test_name"
        timeout "$timeout_duration" bash "$test_path" 2>&1 | tee -a "$test_log"
        local test_result=${PIPESTATUS[0]}
    else
        log_info "Running test: $test_name"
        timeout "$timeout_duration" bash "$test_path" >> "$test_log" 2>&1
        local test_result=$?
    fi
    
    local test_end_time=$(date +%s)
    local test_duration=$((test_end_time - test_start_time))
    
    # Log test results
    echo "=== TEST: $test_name ===" >> "$TEST_SUITE_FILE"
    echo "Duration: ${test_duration}s" >> "$TEST_SUITE_FILE"
    echo "Result: $([ $test_result -eq 0 ] && echo "PASS" || echo "FAIL")" >> "$TEST_SUITE_FILE"
    echo "Log: $test_log" >> "$TEST_SUITE_FILE"
    echo "" >> "$TEST_SUITE_FILE"
    
    log_info "Test $test_name completed in ${test_duration}s"
    
    return $test_result
}

# Run post-test analysis
run_post_test_analysis() {
    local test_suite="$1"
    
    log_info "Running post-test analysis..."
    
    # Analyze test results
    local total_assertions=$(grep -c "PASS\|FAIL" "$TEST_RESULTS_DIR"/*.log 2>/dev/null || echo "0")
    local passed_assertions=$(grep -c "PASS" "$TEST_RESULTS_DIR"/*.log 2>/dev/null || echo "0")
    local failed_assertions=$(grep -c "FAIL" "$TEST_RESULTS_DIR"/*.log 2>/dev/null || echo "0")
    
    echo "=== POST-TEST ANALYSIS ===" >> "$TEST_SUITE_FILE"
    echo "Total assertions: $total_assertions" >> "$TEST_SUITE_FILE"
    echo "Passed assertions: $passed_assertions" >> "$TEST_SUITE_FILE"
    echo "Failed assertions: $failed_assertions" >> "$TEST_SUITE_FILE"
    echo "" >> "$TEST_SUITE_FILE"
    
    # Performance analysis
    if [ "$test_suite" = "performance" ] || [ "$test_suite" = "full" ]; then
        run_performance_analysis
    fi
    
    # Generate recommendations
    generate_test_recommendations
    
    log_success "Post-test analysis completed"
}

# Run performance analysis
run_performance_analysis() {
    log_info "Running performance analysis..."
    
    # Analyze test execution times
    echo "=== PERFORMANCE ANALYSIS ===" >> "$TEST_SUITE_FILE"
    
    # Get test durations from logs
    grep "Duration:" "$TEST_SUITE_FILE" | while read -r line; do
        echo "$line" >> "$TEST_SUITE_FILE"
    done
    
    # System resource usage
    echo "Memory usage:" >> "$TEST_SUITE_FILE"
    free -h >> "$TEST_SUITE_FILE" 2>/dev/null || echo "Memory info unavailable" >> "$TEST_SUITE_FILE"
    
    echo "Disk usage:" >> "$TEST_SUITE_FILE"
    df -h /workspace >> "$TEST_SUITE_FILE" 2>/dev/null || echo "Disk info unavailable" >> "$TEST_SUITE_FILE"
    
    echo "" >> "$TEST_SUITE_FILE"
}

# Generate test recommendations
generate_test_recommendations() {
    log_info "Generating test recommendations..."
    
    echo "=== RECOMMENDATIONS ===" >> "$TEST_SUITE_FILE"
    
    # Analyze failed tests for patterns
    local failed_pattern_count=$(grep -c "FAIL" "$TEST_SUITE_FILE" 2>/dev/null || echo "0")
    
    if [ "$failed_pattern_count" -gt 0 ]; then
        echo "- Review failed tests and fix underlying issues" >> "$TEST_SUITE_FILE"
        echo "- Consider increasing test timeouts if tests are timing out" >> "$TEST_SUITE_FILE"
    fi
    
    # Check for performance issues
    local slow_tests=$(grep "Duration:" "$TEST_SUITE_FILE" | awk '$2 > 300 {print $0}' | wc -l)
    if [ "$slow_tests" -gt 0 ]; then
        echo "- Optimize slow-running tests (>5 minutes)" >> "$TEST_SUITE_FILE"
    fi
    
    # Memory recommendations
    local memory_gb=$(free -g | awk '/^Mem:/{print $2}')
    if [ "$memory_gb" -lt 4 ]; then
        echo "- Consider increasing available memory (current: ${memory_gb}GB, recommended: 4GB+)" >> "$TEST_SUITE_FILE"
    fi
    
    echo "" >> "$TEST_SUITE_FILE"
}

# Generate final report
generate_final_report() {
    local total_tests="$1"
    local passed_tests="$2"
    local failed_tests="$3"
    
    log_info "Generating final test report..."
    
    echo "=== FINAL REPORT ===" >> "$TEST_SUITE_FILE"
    echo "Total test files: $total_tests" >> "$TEST_SUITE_FILE"
    echo "Passed: $passed_tests" >> "$TEST_SUITE_FILE"
    echo "Failed: $failed_tests" >> "$TEST_SUITE_FILE"
    echo "Success rate: $(( passed_tests * 100 / total_tests ))%" >> "$TEST_SUITE_FILE"
    echo "Completed: $(date)" >> "$TEST_SUITE_FILE"
    echo "" >> "$TEST_SUITE_FILE"
    
    # CI mode specific output
    if [ "$CI_MODE" = "true" ]; then
        echo "::notice::Test suite completed with $passed_tests/$total_tests tests passing"
        if [ $failed_tests -gt 0 ]; then
            echo "::error::$failed_tests tests failed"
        fi
    fi
    
    log_info "Test report saved to: $TEST_SUITE_FILE"
}

# Show usage information
usage() {
    cat << EOF
Usage: $0 [SUITE] [VERBOSE] [FAIL_FAST]

Run automated test suites for the dev container.

Arguments:
  SUITE       Test suite to run (default: quick)
              Available: ${!TEST_SUITES[@]}
  
  VERBOSE     Enable verbose output (default: false)
              Values: true, false
  
  FAIL_FAST   Stop on first test failure (default: false)
              Values: true, false

Environment Variables:
  CI=true                Enable CI mode
  TEST_MODE=automated    Set test mode

Test Suites:
  quick       Fast environment validation only
  installation Validate package installations
  full        Complete test suite
  performance Full suite with performance analysis
  ci          CI-optimized test suite

Examples:
  $0                     # Run quick test suite
  $0 full               # Run full test suite
  $0 installation true  # Run installation tests with verbose output
  $0 ci false true     # Run CI tests with fail-fast mode

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