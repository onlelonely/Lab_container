#!/bin/bash

# Configuration management system for dev container
set -euo pipefail

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UTILS_DIR="$SCRIPT_DIR/../utils"
CONFIG_DIR="$SCRIPT_DIR/../../configs"

# Source utility functions
source "$UTILS_DIR/logging.sh"
source "$UTILS_DIR/validation.sh"

# Script name for logging
readonly SCRIPT_NAME="config-manager.sh"

# Configuration paths
readonly CONFIG_BACKUP_DIR="/workspace/logs/config-backups"
readonly CONFIG_TEMPLATES_DIR="$CONFIG_DIR/templates"
readonly USER_CONFIG_DIR="/workspace/.devcontainer-config"

# Create necessary directories
mkdir -p "$CONFIG_BACKUP_DIR" "$USER_CONFIG_DIR"

# Main configuration management function
main() {
    log_script_start "$SCRIPT_NAME"
    
    # Parse command line arguments
    local action="${1:-help}"
    local target="${2:-all}"
    local value="${3:-}"
    
    case "$action" in
        "list")
            list_configurations "$target"
            ;;
        "get")
            get_configuration "$target"
            ;;
        "set")
            set_configuration "$target" "$value"
            ;;
        "reset")
            reset_configuration "$target"
            ;;
        "backup")
            backup_configurations "$target"
            ;;
        "restore")
            restore_configurations "$target"
            ;;
        "validate")
            validate_configurations "$target"
            ;;
        "template")
            create_configuration_template "$target"
            ;;
        "export")
            export_configurations "$target"
            ;;
        "import")
            import_configurations "$target"
            ;;
        "help")
            show_usage
            ;;
        *)
            log_error "Unknown action: $action"
            show_usage
            exit 1
            ;;
    esac
    
    log_script_end "$SCRIPT_NAME" 0
}

# List configurations
list_configurations() {
    local target="$1"
    
    log_info "Listing configurations for: $target"
    
    case "$target" in
        "all"|"environment")
            log_info "Environment variables:"
            list_environment_variables
            ;;
    esac
    
    case "$target" in
        "all"|"python")
            log_info "Python packages:"
            list_python_packages
            ;;
    esac
    
    case "$target" in
        "all"|"r")
            log_info "R packages:"
            list_r_packages
            ;;
    esac
    
    case "$target" in
        "all"|"conda")
            log_info "Conda packages:"
            list_conda_packages
            ;;
    esac
}

# List environment variables
list_environment_variables() {
    local env_vars=(
        "PYTHON_VERSION"
        "R_VERSION"
        "WORKSPACE_DIR"
        "INSTALL_BIOCONDUCTOR"
        "INSTALL_EXTENDED_PACKAGES"
        "PARALLEL_JOBS"
        "ENABLE_LOGGING"
        "LOG_LEVEL"
    )
    
    for var in "${env_vars[@]}"; do
        local value="${!var:-<not set>}"
        echo "  $var = $value"
    done
}

# List Python packages
list_python_packages() {
    local config_files=(
        "requirements-core.txt"
        "requirements-bioinformatics.txt"
        "requirements-extended.txt"
    )
    
    for file in "${config_files[@]}"; do
        local file_path="$CONFIG_DIR/python/$file"
        if [ -f "$file_path" ]; then
            echo "  $file:"
            grep -v "^#" "$file_path" | grep -v "^$" | sed 's/^/    /'
        fi
    done
}

# List R packages
list_r_packages() {
    local config_files=(
        "core-packages.R"
        "bioconductor-packages.R"
        "extended-packages.R"
    )
    
    for file in "${config_files[@]}"; do
        local file_path="$CONFIG_DIR/r/$file"
        if [ -f "$file_path" ]; then
            echo "  $file:"
            grep -o '"[^"]*"' "$file_path" | sed 's/"//g' | sed 's/^/    /'
        fi
    done
}

# List Conda packages
list_conda_packages() {
    local config_files=(
        "environment-core.yml"
        "environment-extended.yml"
    )
    
    for file in "${config_files[@]}"; do
        local file_path="$CONFIG_DIR/conda/$file"
        if [ -f "$file_path" ]; then
            echo "  $file:"
            if command -v yq > /dev/null 2>&1; then
                yq e '.dependencies[]' "$file_path" | sed 's/^/    /'
            else
                grep -A 20 "dependencies:" "$file_path" | grep "  -" | sed 's/^/    /'
            fi
        fi
    done
}

# Get configuration value
get_configuration() {
    local target="$1"
    
    case "$target" in
        "python-version")
            echo "${PYTHON_VERSION:-3.11}"
            ;;
        "r-version")
            echo "${R_VERSION:-4.3}"
            ;;
        "parallel-jobs")
            echo "${PARALLEL_JOBS:-4}"
            ;;
        "install-bioconductor")
            echo "${INSTALL_BIOCONDUCTOR:-true}"
            ;;
        "install-extended")
            echo "${INSTALL_EXTENDED_PACKAGES:-false}"
            ;;
        "log-level")
            echo "${LOG_LEVEL:-INFO}"
            ;;
        *)
            log_error "Unknown configuration: $target"
            return 1
            ;;
    esac
}

# Set configuration value
set_configuration() {
    local target="$1"
    local value="$2"
    
    if [ -z "$value" ]; then
        log_error "Value required for setting configuration"
        return 1
    fi
    
    log_info "Setting $target = $value"
    
    # Create user configuration file if it doesn't exist
    local user_config="$USER_CONFIG_DIR/environment.conf"
    touch "$user_config"
    
    case "$target" in
        "python-version")
            set_env_var "PYTHON_VERSION" "$value"
            ;;
        "r-version")
            set_env_var "R_VERSION" "$value"
            ;;
        "parallel-jobs")
            if [[ "$value" =~ ^[0-9]+$ ]] && [ "$value" -gt 0 ] && [ "$value" -le 32 ]; then
                set_env_var "PARALLEL_JOBS" "$value"
            else
                log_error "Invalid parallel jobs value: $value (must be 1-32)"
                return 1
            fi
            ;;
        "install-bioconductor")
            if [[ "$value" =~ ^(true|false)$ ]]; then
                set_env_var "INSTALL_BIOCONDUCTOR" "$value"
            else
                log_error "Invalid value: $value (must be true or false)"
                return 1
            fi
            ;;
        "install-extended")
            if [[ "$value" =~ ^(true|false)$ ]]; then
                set_env_var "INSTALL_EXTENDED_PACKAGES" "$value"
            else
                log_error "Invalid value: $value (must be true or false)"
                return 1
            fi
            ;;
        "log-level")
            if [[ "$value" =~ ^(DEBUG|INFO|WARNING|ERROR)$ ]]; then
                set_env_var "LOG_LEVEL" "$value"
            else
                log_error "Invalid log level: $value (must be DEBUG, INFO, WARNING, or ERROR)"
                return 1
            fi
            ;;
        *)
            log_error "Unknown configuration: $target"
            return 1
            ;;
    esac
    
    log_success "Configuration updated: $target = $value"
}

# Set environment variable in user config
set_env_var() {
    local var_name="$1"
    local var_value="$2"
    local user_config="$USER_CONFIG_DIR/environment.conf"
    
    # Remove existing setting
    sed -i "/^$var_name=/d" "$user_config" 2>/dev/null || true
    
    # Add new setting
    echo "$var_name=$var_value" >> "$user_config"
    
    # Export for current session
    export "$var_name=$var_value"
}

# Reset configuration to defaults
reset_configuration() {
    local target="$1"
    
    log_info "Resetting configuration: $target"
    
    case "$target" in
        "all")
            rm -f "$USER_CONFIG_DIR/environment.conf"
            log_success "All configurations reset to defaults"
            ;;
        *)
            local user_config="$USER_CONFIG_DIR/environment.conf"
            if [ -f "$user_config" ]; then
                sed -i "/^$target=/d" "$user_config" 2>/dev/null || true
                log_success "Configuration reset: $target"
            else
                log_info "No user configuration to reset"
            fi
            ;;
    esac
}

# Backup configurations
backup_configurations() {
    local target="$1"
    local backup_name="backup-$(date +%Y%m%d-%H%M%S)"
    local backup_path="$CONFIG_BACKUP_DIR/$backup_name"
    
    log_info "Creating configuration backup: $backup_name"
    
    mkdir -p "$backup_path"
    
    case "$target" in
        "all")
            cp -r "$CONFIG_DIR"/* "$backup_path/"
            if [ -f "$USER_CONFIG_DIR/environment.conf" ]; then
                cp "$USER_CONFIG_DIR/environment.conf" "$backup_path/"
            fi
            ;;
        "python")
            cp -r "$CONFIG_DIR/python" "$backup_path/"
            ;;
        "r")
            cp -r "$CONFIG_DIR/r" "$backup_path/"
            ;;
        "conda")
            cp -r "$CONFIG_DIR/conda" "$backup_path/"
            ;;
        *)
            log_error "Unknown backup target: $target"
            return 1
            ;;
    esac
    
    log_success "Backup created: $backup_path"
}

# Restore configurations
restore_configurations() {
    local backup_name="$1"
    
    if [ -z "$backup_name" ]; then
        log_info "Available backups:"
        ls -1 "$CONFIG_BACKUP_DIR" 2>/dev/null || log_info "No backups found"
        return 0
    fi
    
    local backup_path="$CONFIG_BACKUP_DIR/$backup_name"
    
    if [ ! -d "$backup_path" ]; then
        log_error "Backup not found: $backup_name"
        return 1
    fi
    
    log_info "Restoring configuration from backup: $backup_name"
    
    # Create current backup before restoring
    backup_configurations "all"
    
    # Restore from backup
    cp -r "$backup_path"/* "$CONFIG_DIR/"
    
    if [ -f "$backup_path/environment.conf" ]; then
        cp "$backup_path/environment.conf" "$USER_CONFIG_DIR/"
    fi
    
    log_success "Configuration restored from: $backup_name"
}

# Validate configurations
validate_configurations() {
    local target="$1"
    
    log_info "Validating configurations: $target"
    
    local validation_failed=0
    
    case "$target" in
        "all"|"python")
            if ! validate_python_configs; then
                validation_failed=1
            fi
            ;;
    esac
    
    case "$target" in
        "all"|"r")
            if ! validate_r_configs; then
                validation_failed=1
            fi
            ;;
    esac
    
    case "$target" in
        "all"|"conda")
            if ! validate_conda_configs; then
                validation_failed=1
            fi
            ;;
    esac
    
    if [ $validation_failed -eq 0 ]; then
        log_success "All configurations are valid"
        return 0
    else
        log_error "Configuration validation failed"
        return 1
    fi
}

# Validate Python configurations
validate_python_configs() {
    local config_files=(
        "requirements-core.txt"
        "requirements-bioinformatics.txt"
        "requirements-extended.txt"
    )
    
    local validation_passed=0
    
    for file in "${config_files[@]}"; do
        local file_path="$CONFIG_DIR/python/$file"
        if [ -f "$file_path" ]; then
            # Check for valid package format
            if grep -q "^[a-zA-Z0-9_-]" "$file_path"; then
                log_success "✓ $file is valid"
            else
                log_error "✗ $file has invalid format"
                validation_passed=1
            fi
        else
            log_error "✗ $file not found"
            validation_passed=1
        fi
    done
    
    return $validation_passed
}

# Validate R configurations
validate_r_configs() {
    local config_files=(
        "core-packages.R"
        "bioconductor-packages.R"
        "extended-packages.R"
    )
    
    local validation_passed=0
    
    for file in "${config_files[@]}"; do
        local file_path="$CONFIG_DIR/r/$file"
        if [ -f "$file_path" ]; then
            # Check for valid R syntax
            if R --slave -e "source('$file_path')" 2>/dev/null; then
                log_success "✓ $file is valid"
            else
                log_error "✗ $file has invalid R syntax"
                validation_passed=1
            fi
        else
            log_error "✗ $file not found"
            validation_passed=1
        fi
    done
    
    return $validation_passed
}

# Validate Conda configurations
validate_conda_configs() {
    local config_files=(
        "environment-core.yml"
        "environment-extended.yml"
    )
    
    local validation_passed=0
    
    for file in "${config_files[@]}"; do
        local file_path="$CONFIG_DIR/conda/$file"
        if [ -f "$file_path" ]; then
            # Check YAML syntax
            if command -v python3 > /dev/null 2>&1; then
                if python3 -c "import yaml; yaml.safe_load(open('$file_path'))" 2>/dev/null; then
                    log_success "✓ $file is valid YAML"
                else
                    log_error "✗ $file has invalid YAML syntax"
                    validation_passed=1
                fi
            else
                log_warning "Python3 not available for YAML validation"
            fi
        else
            log_error "✗ $file not found"
            validation_passed=1
        fi
    done
    
    return $validation_passed
}

# Create configuration template
create_configuration_template() {
    local target="$1"
    
    log_info "Creating configuration template: $target"
    
    mkdir -p "$CONFIG_TEMPLATES_DIR"
    
    case "$target" in
        "environment")
            create_environment_template
            ;;
        "python")
            create_python_template
            ;;
        "r")
            create_r_template
            ;;
        "conda")
            create_conda_template
            ;;
        *)
            log_error "Unknown template target: $target"
            return 1
            ;;
    esac
}

# Create environment template
create_environment_template() {
    local template_file="$CONFIG_TEMPLATES_DIR/environment.conf.template"
    
    cat > "$template_file" << 'EOF'
# Dev Container Environment Configuration Template
# Copy this file to ~/.devcontainer-config/environment.conf and customize

# Core Configuration
PYTHON_VERSION=3.11
R_VERSION=4.3
WORKSPACE_DIR=/workspace

# Package Installation Settings
INSTALL_BIOCONDUCTOR=true
INSTALL_EXTENDED_PACKAGES=false
PARALLEL_JOBS=4

# Logging Configuration
ENABLE_LOGGING=true
LOG_LEVEL=INFO

# Development Settings
ENABLE_VALIDATION=true
EOF
    
    log_success "Environment template created: $template_file"
}

# Export configurations
export_configurations() {
    local target="$1"
    local export_file="/workspace/devcontainer-config-export.tar.gz"
    
    log_info "Exporting configurations to: $export_file"
    
    local temp_dir=$(mktemp -d)
    
    case "$target" in
        "all")
            cp -r "$CONFIG_DIR" "$temp_dir/"
            if [ -f "$USER_CONFIG_DIR/environment.conf" ]; then
                cp "$USER_CONFIG_DIR/environment.conf" "$temp_dir/"
            fi
            ;;
        *)
            log_error "Export target not supported: $target"
            return 1
            ;;
    esac
    
    tar -czf "$export_file" -C "$temp_dir" .
    rm -rf "$temp_dir"
    
    log_success "Configurations exported to: $export_file"
}

# Import configurations
import_configurations() {
    local import_file="$1"
    
    if [ -z "$import_file" ] || [ ! -f "$import_file" ]; then
        log_error "Import file required and must exist"
        return 1
    fi
    
    log_info "Importing configurations from: $import_file"
    
    # Create backup before importing
    backup_configurations "all"
    
    # Extract to temporary directory
    local temp_dir=$(mktemp -d)
    tar -xzf "$import_file" -C "$temp_dir"
    
    # Import configurations
    if [ -d "$temp_dir/configs" ]; then
        cp -r "$temp_dir/configs"/* "$CONFIG_DIR/"
    fi
    
    if [ -f "$temp_dir/environment.conf" ]; then
        cp "$temp_dir/environment.conf" "$USER_CONFIG_DIR/"
    fi
    
    rm -rf "$temp_dir"
    
    log_success "Configurations imported successfully"
}

# Show usage information
show_usage() {
    cat << EOF
Usage: $0 ACTION [TARGET] [VALUE]

Configuration management system for dev container.

Actions:
  list [TARGET]           List configurations
  get TARGET              Get configuration value
  set TARGET VALUE        Set configuration value
  reset [TARGET]          Reset configuration to defaults
  backup [TARGET]         Backup configurations
  restore [BACKUP_NAME]   Restore from backup
  validate [TARGET]       Validate configurations
  template TARGET         Create configuration template
  export TARGET           Export configurations
  import FILE             Import configurations

Targets:
  all                     All configurations
  environment             Environment variables
  python                  Python packages
  r                       R packages
  conda                   Conda packages

Configuration Keys:
  python-version          Python version (e.g., 3.11)
  r-version              R version (e.g., 4.3)
  parallel-jobs          Number of parallel jobs (1-32)
  install-bioconductor   Install Bioconductor (true/false)
  install-extended       Install extended packages (true/false)
  log-level              Log level (DEBUG/INFO/WARNING/ERROR)

Examples:
  $0 list                        # List all configurations
  $0 get python-version         # Get Python version
  $0 set parallel-jobs 8        # Set parallel jobs to 8
  $0 backup all                 # Backup all configurations
  $0 validate python            # Validate Python configs

EOF
}

# Error handling
trap 'log_script_end "$SCRIPT_NAME" $?' EXIT

# Run main function
main "$@"