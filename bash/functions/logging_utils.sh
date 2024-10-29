#!/bin/bash

# Load required settings for logging_utils
source "$HOME/lab_utils/bash/config/project_config.sh"
# Advanced Logging Functions for Bash Scripts
#
# Script: 002_logging_functions.sh
# Description: A set of functions for consistent logging across Bash scripts
# Author: Your Name
# Date: 2024-10-18

# Function: get_script_name
# Purpose: Extract the full path of the current script
# Parameters: None
# Return: Script path
get_script_name() {
  echo "$(readlink -f "${BASH_SOURCE[0]}")"
}

# Function: get_script_dir
# Purpose: Extract the directory of the current script
# Parameters: None
# Return: Script directory
get_script_dir() {
  echo "$(dirname "$(get_script_name)")"
}

# Function: get_script_basename
# Purpose: Extract the basename of the current script without extension
# Parameters: None
# Return: Script basename
get_script_basename() {
  echo "$(basename "$(get_script_name)" .sh)"
}

#' Initialize Logging System
#' @param script_name Character Name of the calling script
#' @param log_dir Character Optional custom log directory
#' @return String Path to log file
initialize_logging() {
    local script_name="${1:-$(basename "${BASH_SOURCE[1]}" )}"
    local log_dir="${2:-${PROJECT_CONFIG[DEFAULT_LOG_ROOT]}}"
    local -n return_var=$3
    log_info "Logging initialized"
    
    # Ensure log directory exists
    mkdir -p "$log_dir/$(date +%Y-%m)"
    
    # Generate log file path
    return_var="$log_dir/$(date +%Y-%m)/$(date +%Y-%m-%d)_${script_name}.log"
    #
    # Add run separator and count runs
    if [[ -f "$return_var" ]]; then
        # Count existing runs
        local run_count=$(grep -c "=== New Run ===" "$return_var")
        run_count=$((run_count + 1))
        # Add separator with run count
        echo -e "\n=== New Run === (#$run_count) === $(date +'%Y-%m-%d %H:%M:%S') ===\n" >> "$return_var"
    else
        # First run - create file and add initial header
        echo "=== New Run === (#1) === $(date +'%Y-%m-%d %H:%M:%S') ===" > "$return_var"
    fi
    # Initialize log file with headers
    log_system_info "$return_var"
    log_git_info "$return_var"
}

#' Log System Information
#' @param log_file Character Path to log file
#' @return None
log_system_info() {
    local log_file="$1"
    
    log_message "INFO" "System Information:" "$log_file"
    log_message "INFO" "  Bash: $BASH_VERSION" "$log_file"
    log_message "INFO" "  Host: $(hostname)" "$log_file"
    log_message "INFO" "  User: $USER" "$log_file"
    log_message "INFO" "  PWD:  $PWD" "$log_file"
}

#' Log Git Information
#' @param log_file Character Path to log file
#' @return None
log_git_info() {
    local log_file="$1"
    
    if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
        log_message "INFO" "Git Information:" "$log_file"
        log_message "INFO" "  Branch: $(git rev-parse --abbrev-ref HEAD)" "$log_file"
        log_message "INFO" "  Commit: $(git rev-parse HEAD)" "$log_file"
    else
        log_message "WARNING" "Not in a git repository" "$log_file"
    fi
}

#' Log Message
#' @param level Character Log level
#' @param message Character Message to log
#' @param log_file Character Path to log file
#' @return None
log_message() {
    local level="$1"
    local message="$2"
    local log_file="$3"
    
    # Validate log level
    if [[ ! " ${PROJECT_CONFIG[LOG_LEVELS]} " =~ " ${level} " ]]; then
        echo "Invalid log level: $level" >&2
        return 1
    fi
    
    # Format message
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    local log_entry="[${timestamp}] [${level}] ${message}"
    
    # Output to console
    echo "$log_entry"
    
    # Write to log file if provided
    if [[ -n "$log_file" ]]; then
        echo "$log_entry" >> "$log_file"
    fi
}

#' Convenience Logging Functions
for level in ${PROJECT_CONFIG[LOG_LEVELS]}; do
    eval "log_${level,,}() { log_message \"$level\" \"\$1\" \"\$2\"; }"
done