#!/bin/bash

# Unified logging system for dev container
set -euo pipefail

# Prevent multiple sourcing issues with readonly variables
if [ -z "${LOG_DIR:-}" ]; then
    readonly LOG_DIR="/workspace/logs"
fi

if [ -z "${LOG_FILE:-}" ]; then
    readonly LOG_FILE="${LOG_DIR}/devcontainer-$(date +%Y%m%d).log"
fi

# Ensure log directory exists
mkdir -p "$LOG_DIR"

# Log levels (with existence checks)
if [ -z "${LOG_LEVEL_INFO:-}" ]; then
    readonly LOG_LEVEL_INFO="INFO"
fi

if [ -z "${LOG_LEVEL_WARNING:-}" ]; then
    readonly LOG_LEVEL_WARNING="WARNING"
fi

if [ -z "${LOG_LEVEL_ERROR:-}" ]; then
    readonly LOG_LEVEL_ERROR="ERROR"
fi

if [ -z "${LOG_LEVEL_SUCCESS:-}" ]; then
    readonly LOG_LEVEL_SUCCESS="SUCCESS"
fi

# Color codes for terminal output (with existence checks)
if [ -z "${COLOR_RESET:-}" ]; then
    readonly COLOR_RESET='\033[0m'
fi

if [ -z "${COLOR_INFO:-}" ]; then
    readonly COLOR_INFO='\033[0;34m'      # Blue
fi

if [ -z "${COLOR_WARNING:-}" ]; then
    readonly COLOR_WARNING='\033[0;33m'   # Yellow
fi

if [ -z "${COLOR_ERROR:-}" ]; then
    readonly COLOR_ERROR='\033[0;31m'     # Red
fi

if [ -z "${COLOR_SUCCESS:-}" ]; then
    readonly COLOR_SUCCESS='\033[0;32m'   # Green
fi

# Log function with level
log_with_level() {
    local level="$1"
    local message="$2"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    local log_entry="[$timestamp] [$level] $message"
    
    # Write to log file
    echo "$log_entry" >> "$LOG_FILE"
    
    # Output to terminal with color
    case "$level" in
        "$LOG_LEVEL_INFO")
            echo -e "${COLOR_INFO}[INFO]${COLOR_RESET} $message" >&2
            ;;
        "$LOG_LEVEL_WARNING")
            echo -e "${COLOR_WARNING}[WARNING]${COLOR_RESET} $message" >&2
            ;;
        "$LOG_LEVEL_ERROR")
            echo -e "${COLOR_ERROR}[ERROR]${COLOR_RESET} $message" >&2
            ;;
        "$LOG_LEVEL_SUCCESS")
            echo -e "${COLOR_SUCCESS}[SUCCESS]${COLOR_RESET} $message" >&2
            ;;
        *)
            echo "[$level] $message" >&2
            ;;
    esac
}

# Convenience functions
log_info() {
    log_with_level "$LOG_LEVEL_INFO" "$1"
}

log_warning() {
    log_with_level "$LOG_LEVEL_WARNING" "$1"
}

log_error() {
    log_with_level "$LOG_LEVEL_ERROR" "$1"
}

log_success() {
    log_with_level "$LOG_LEVEL_SUCCESS" "$1"
}

# Function to log command execution
log_command() {
    local command="$1"
    log_info "Executing: $command"
}

# Function to log script start
log_script_start() {
    local script_name="$1"
    log_info "Starting script: $script_name"
    log_info "Log file: $LOG_FILE"
}

# Function to log script end
log_script_end() {
    local script_name="$1"
    local exit_code="$2"
    if [ "$exit_code" -eq 0 ]; then
        log_success "Script completed successfully: $script_name"
    else
        log_error "Script failed with exit code $exit_code: $script_name"
    fi
}

# Export functions for use in other scripts
export -f log_with_level log_info log_warning log_error log_success log_command log_script_start log_script_end