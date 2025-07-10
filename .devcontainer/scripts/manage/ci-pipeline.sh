#!/bin/bash

# CI/CD pipeline automation for dev container
set -euo pipefail

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UTILS_DIR="$SCRIPT_DIR/../utils"

# Source utility functions
source "$UTILS_DIR/logging.sh"
source "$UTILS_DIR/validation.sh"

# Script name for logging
readonly SCRIPT_NAME="ci-pipeline.sh"

# CI configuration
readonly CI_RESULTS_DIR="/workspace/logs/ci-results"
readonly CI_REPORT_FILE="$CI_RESULTS_DIR/ci-pipeline-$(date +%Y%m%d-%H%M%S).log"

# Create CI results directory
mkdir -p "$CI_RESULTS_DIR"

# Main CI pipeline function
main() {
    log_script_start "$SCRIPT_NAME"
    
    # Parse command line arguments
    local stage="${1:-all}"
    local environment="${2:-development}"
    
    log_info "Running CI pipeline stage: $stage"
    log_info "Environment: $environment"
    
    # Initialize CI report
    echo "=== CI/CD PIPELINE REPORT ===" > "$CI_REPORT_FILE"
    echo "Stage: $stage" >> "$CI_REPORT_FILE"
    echo "Environment: $environment" >> "$CI_REPORT_FILE"
    echo "Started: $(date)" >> "$CI_REPORT_FILE"
    echo "" >> "$CI_REPORT_FILE"
    
    # Set CI environment
    export CI=true
    export CI_ENVIRONMENT="$environment"
    
    # Run pipeline stages
    case "$stage" in
        "validate")
            run_validation_stage
            ;;
        "build")
            run_build_stage
            ;;
        "test")
            run_test_stage
            ;;
        "deploy")
            run_deploy_stage "$environment"
            ;;
        "all")
            run_validation_stage && \
            run_build_stage && \
            run_test_stage && \
            run_deploy_stage "$environment"
            ;;
        *)
            log_error "Unknown CI stage: $stage"
            log_info "Available stages: validate, build, test, deploy, all"
            exit 1
            ;;
    esac
    
    # Generate CI summary
    generate_ci_summary
    
    log_script_end "$SCRIPT_NAME" 0
}

# Validation stage
run_validation_stage() {
    log_info "=== CI Stage: Validation ==="
    
    echo "=== VALIDATION STAGE ===" >> "$CI_REPORT_FILE"
    local stage_start_time=$(date +%s)
    
    # Validate configuration files
    log_info "Validating configuration files..."
    local config_validation=0
    
    # Check Dockerfile
    if [ -f "$SCRIPT_DIR/../../Dockerfile" ]; then
        log_success "✓ Dockerfile found"
        echo "PASS: Dockerfile exists" >> "$CI_REPORT_FILE"
    else
        log_error "✗ Dockerfile not found"
        echo "FAIL: Dockerfile missing" >> "$CI_REPORT_FILE"
        config_validation=1
    fi
    
    # Check docker-compose.yml
    if [ -f "$SCRIPT_DIR/../../docker-compose.yml" ]; then
        log_success "✓ docker-compose.yml found"
        echo "PASS: docker-compose.yml exists" >> "$CI_REPORT_FILE"
    else
        log_error "✗ docker-compose.yml not found"
        echo "FAIL: docker-compose.yml missing" >> "$CI_REPORT_FILE"
        config_validation=1
    fi
    
    # Check devcontainer.json
    if [ -f "$SCRIPT_DIR/../../devcontainer.json" ]; then
        log_success "✓ devcontainer.json found"
        echo "PASS: devcontainer.json exists" >> "$CI_REPORT_FILE"
        
        # Validate JSON syntax
        if command -v jq > /dev/null 2>&1; then
            if jq empty "$SCRIPT_DIR/../../devcontainer.json" 2>/dev/null; then
                log_success "✓ devcontainer.json syntax valid"
                echo "PASS: devcontainer.json syntax valid" >> "$CI_REPORT_FILE"
            else
                log_error "✗ devcontainer.json syntax invalid"
                echo "FAIL: devcontainer.json syntax invalid" >> "$CI_REPORT_FILE"
                config_validation=1
            fi
        fi
    else
        log_error "✗ devcontainer.json not found"
        echo "FAIL: devcontainer.json missing" >> "$CI_REPORT_FILE"
        config_validation=1
    fi
    
    # Validate configuration directories
    local config_dirs=("configs/python" "configs/r" "configs/conda" "scripts/utils" "tests")
    for dir in "${config_dirs[@]}"; do
        if [ -d "$SCRIPT_DIR/../../$dir" ]; then
            log_success "✓ Directory found: $dir"
            echo "PASS: Directory exists: $dir" >> "$CI_REPORT_FILE"
        else
            log_error "✗ Directory missing: $dir"
            echo "FAIL: Directory missing: $dir" >> "$CI_REPORT_FILE"
            config_validation=1
        fi
    done
    
    local stage_end_time=$(date +%s)
    local stage_duration=$((stage_end_time - stage_start_time))
    
    echo "Duration: ${stage_duration}s" >> "$CI_REPORT_FILE"
    echo "Result: $([ $config_validation -eq 0 ] && echo "PASS" || echo "FAIL")" >> "$CI_REPORT_FILE"
    echo "" >> "$CI_REPORT_FILE"
    
    if [ $config_validation -eq 0 ]; then
        log_success "Validation stage passed"
        return 0
    else
        log_error "Validation stage failed"
        return 1
    fi
}

# Build stage
run_build_stage() {
    log_info "=== CI Stage: Build ==="
    
    echo "=== BUILD STAGE ===" >> "$CI_REPORT_FILE"
    local stage_start_time=$(date +%s)
    
    # Simulate build process (in real CI, this would build the Docker image)
    log_info "Simulating Docker image build..."
    
    # Check if Docker is available
    if command -v docker > /dev/null 2>&1; then
        log_info "Docker available, could build image"
        echo "INFO: Docker available for building" >> "$CI_REPORT_FILE"
        
        # Check if we can access Docker daemon
        if docker info > /dev/null 2>&1; then
            log_info "Docker daemon accessible"
            echo "INFO: Docker daemon accessible" >> "$CI_REPORT_FILE"
            
            # In a real CI environment, you would run:
            # docker build -t devcontainer -f .devcontainer/Dockerfile .
            log_info "Build simulation: docker build would run here"
            echo "SIMULATION: Docker build process" >> "$CI_REPORT_FILE"
        else
            log_warning "Docker daemon not accessible (expected in dev container)"
            echo "WARNING: Docker daemon not accessible" >> "$CI_REPORT_FILE"
        fi
    else
        log_warning "Docker not available (expected in dev container)"
        echo "WARNING: Docker not available" >> "$CI_REPORT_FILE"
    fi
    
    # Validate build artifacts
    log_info "Validating build artifacts..."
    local required_scripts=(
        "scripts/install/install-core.sh"
        "scripts/install/install-extended.sh"
        "scripts/install/install-parallel.sh"
        "scripts/utils/logging.sh"
        "scripts/utils/retry-logic.sh"
        "scripts/utils/validation.sh"
    )
    
    local build_validation=0
    for script in "${required_scripts[@]}"; do
        if [ -f "$SCRIPT_DIR/../../$script" ] && [ -x "$SCRIPT_DIR/../../$script" ]; then
            log_success "✓ Script available and executable: $script"
            echo "PASS: Script ready: $script" >> "$CI_REPORT_FILE"
        else
            log_error "✗ Script missing or not executable: $script"
            echo "FAIL: Script issue: $script" >> "$CI_REPORT_FILE"
            build_validation=1
        fi
    done
    
    local stage_end_time=$(date +%s)
    local stage_duration=$((stage_end_time - stage_start_time))
    
    echo "Duration: ${stage_duration}s" >> "$CI_REPORT_FILE"
    echo "Result: $([ $build_validation -eq 0 ] && echo "PASS" || echo "FAIL")" >> "$CI_REPORT_FILE"
    echo "" >> "$CI_REPORT_FILE"
    
    if [ $build_validation -eq 0 ]; then
        log_success "Build stage passed"
        return 0
    else
        log_error "Build stage failed"
        return 1
    fi
}

# Test stage
run_test_stage() {
    log_info "=== CI Stage: Test ==="
    
    echo "=== TEST STAGE ===" >> "$CI_REPORT_FILE"
    local stage_start_time=$(date +%s)
    
    # Run automated test suite
    log_info "Running automated test suite..."
    
    # Run CI test suite
    if "$SCRIPT_DIR/test-runner.sh" ci false true; then
        log_success "CI test suite passed"
        echo "PASS: CI test suite" >> "$CI_REPORT_FILE"
        local test_result=0
    else
        log_error "CI test suite failed"
        echo "FAIL: CI test suite" >> "$CI_REPORT_FILE"
        local test_result=1
    fi
    
    # Run security checks
    log_info "Running security checks..."
    local security_result=0
    
    # Check for secrets in configuration files
    if grep -r -i "password\|secret\|key\|token" "$SCRIPT_DIR/../../configs/" 2>/dev/null | grep -v -E "(template|example|placeholder)"; then
        log_warning "Potential secrets found in configuration files"
        echo "WARNING: Potential secrets in configs" >> "$CI_REPORT_FILE"
    else
        log_success "✓ No secrets found in configuration files"
        echo "PASS: No secrets in configs" >> "$CI_REPORT_FILE"
    fi
    
    # Check file permissions
    if find "$SCRIPT_DIR/../.." -name "*.sh" ! -perm -u+x 2>/dev/null | grep -q .; then
        log_warning "Some shell scripts are not executable"
        echo "WARNING: Non-executable scripts found" >> "$CI_REPORT_FILE"
    else
        log_success "✓ All shell scripts are executable"
        echo "PASS: Script permissions correct" >> "$CI_REPORT_FILE"
    fi
    
    local stage_end_time=$(date +%s)
    local stage_duration=$((stage_end_time - stage_start_time))
    
    echo "Duration: ${stage_duration}s" >> "$CI_REPORT_FILE"
    echo "Result: $([ $test_result -eq 0 ] && echo "PASS" || echo "FAIL")" >> "$CI_REPORT_FILE"
    echo "" >> "$CI_REPORT_FILE"
    
    if [ $test_result -eq 0 ]; then
        log_success "Test stage passed"
        return 0
    else
        log_error "Test stage failed"
        return 1
    fi
}

# Deploy stage
run_deploy_stage() {
    local environment="$1"
    
    log_info "=== CI Stage: Deploy ($environment) ==="
    
    echo "=== DEPLOY STAGE ($environment) ===" >> "$CI_REPORT_FILE"
    local stage_start_time=$(date +%s)
    
    case "$environment" in
        "development")
            log_info "Development deployment simulation"
            echo "INFO: Development deployment" >> "$CI_REPORT_FILE"
            # In real CI: deploy to development environment
            ;;
        "staging")
            log_info "Staging deployment simulation"
            echo "INFO: Staging deployment" >> "$CI_REPORT_FILE"
            # In real CI: deploy to staging environment
            ;;
        "production")
            log_info "Production deployment simulation"
            echo "INFO: Production deployment" >> "$CI_REPORT_FILE"
            # In real CI: deploy to production environment
            ;;
        *)
            log_error "Unknown environment: $environment"
            echo "FAIL: Unknown environment: $environment" >> "$CI_REPORT_FILE"
            return 1
            ;;
    esac
    
    # Deployment health check
    log_info "Running deployment health check..."
    if validate_all > /dev/null 2>&1; then
        log_success "✓ Deployment health check passed"
        echo "PASS: Deployment health check" >> "$CI_REPORT_FILE"
        local deploy_result=0
    else
        log_error "✗ Deployment health check failed"
        echo "FAIL: Deployment health check" >> "$CI_REPORT_FILE"
        local deploy_result=1
    fi
    
    local stage_end_time=$(date +%s)
    local stage_duration=$((stage_end_time - stage_start_time))
    
    echo "Duration: ${stage_duration}s" >> "$CI_REPORT_FILE"
    echo "Result: $([ $deploy_result -eq 0 ] && echo "PASS" || echo "FAIL")" >> "$CI_REPORT_FILE"
    echo "" >> "$CI_REPORT_FILE"
    
    if [ $deploy_result -eq 0 ]; then
        log_success "Deploy stage passed"
        return 0
    else
        log_error "Deploy stage failed"
        return 1
    fi
}

# Generate CI summary
generate_ci_summary() {
    log_info "Generating CI pipeline summary..."
    
    echo "=== CI PIPELINE SUMMARY ===" >> "$CI_REPORT_FILE"
    
    # Count stage results
    local total_stages=$(grep -c "=== .* STAGE" "$CI_REPORT_FILE")
    local passed_stages=$(grep -c "Result: PASS" "$CI_REPORT_FILE")
    local failed_stages=$(grep -c "Result: FAIL" "$CI_REPORT_FILE")
    
    echo "Total stages: $total_stages" >> "$CI_REPORT_FILE"
    echo "Passed: $passed_stages" >> "$CI_REPORT_FILE"
    echo "Failed: $failed_stages" >> "$CI_REPORT_FILE"
    echo "Success rate: $(( passed_stages * 100 / total_stages ))%" >> "$CI_REPORT_FILE"
    echo "Completed: $(date)" >> "$CI_REPORT_FILE"
    
    log_info "CI pipeline summary:"
    log_info "Stages: $passed_stages/$total_stages passed"
    log_info "Report saved to: $CI_REPORT_FILE"
}

# Show usage information
usage() {
    cat << EOF
Usage: $0 [STAGE] [ENVIRONMENT]

Run CI/CD pipeline stages for the dev container.

Arguments:
  STAGE        Pipeline stage to run (default: all)
               Available: validate, build, test, deploy, all
  
  ENVIRONMENT  Target environment for deployment (default: development)
               Available: development, staging, production

Examples:
  $0                      # Run complete pipeline for development
  $0 test                 # Run test stage only
  $0 deploy staging      # Deploy to staging environment
  $0 all production      # Complete pipeline for production

Pipeline Stages:
  validate    Validate configuration files and structure
  build       Build and validate artifacts
  test        Run automated test suites and security checks
  deploy      Deploy to target environment
  all         Run all stages in sequence

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