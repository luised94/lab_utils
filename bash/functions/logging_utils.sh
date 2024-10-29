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


#' Initialize Logging
#' @param script_name Character Name of the calling script
#' @param log_dir Character Optional log directory
#' @return String Path to log file

#' Initialize Logging
#' @param script_name Character Name of the calling script
#' @param log_dir Character Optional log directory
#' @return String Path to log file
initialize_logging() {
    local script_name="${1:-$(basename "${BASH_SOURCE[1]}" )}"
    local log_dir="${2:-${PROJECT_CONFIG[DEFAULT_LOG_ROOT]}}"
    local log_file
    
    
    echo "DEBUG: Initializing logging with:" >&2
    echo "  script_name: $script_name" >&2
    echo "  log_dir: $log_dir" >&2

    # Ensure log directory exists with verbose error checking
    if ! mkdir -p "$log_dir/$(date +%Y-%m)"; then
        echo "ERROR: Failed to create log directory: $log_dir/$(date +%Y-%m)"
        return 1
    fi
    
    # Generate log file path
    log_file="$log_dir/$(date +%Y-%m)/$(date +%Y-%m-%d)_${script_name}.log"
    echo "[DEBUG] Log file when it is created: ${log_file}" 
    # Ensure log file is writable
    touch "$log_file" 2>/dev/null || {
        echo "ERROR: Cannot create/write to log file: $log_file"
        return 1
    }
    
    # Add run separator and count runs
    if [[ -f "$log_file" ]]; then
        local run_count=$(grep -c "=== New Run ===" "$log_file" || echo "0")
        run_count=$((run_count + 1))
        echo -e "\n=== New Run === (#$run_count) === $(date +'%Y-%m-%d %H:%M:%S') ===\n" >> "$log_file"
    else
        echo "=== New Run === (#1) === $(date +'%Y-%m-%d %H:%M:%S') ===" > "$log_file"
    fi
    
    # Initialize log file with headers (only once)
    log_system_info "$log_file"
    log_git_info "$log_file"
    
    echo -n "$log_file"
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
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    # Validate log level
    if [[ ! " ${PROJECT_CONFIG[LOG_LEVELS]} " =~ " ${level} " ]]; then
        echo "Invalid log level: $level" >&2
        return 1
    fi
    
    # Format message
    local log_entry="[${timestamp}] [${level}] ${message}"
    
    # Output to console
    echo "$log_entry"
    
    # Write to log file if provided and writable
    if [[ -n "$log_file" ]] && [[ -w "$log_file" ]]; then
        echo "$log_entry" >> "$log_file"
    elif [[ -n "$log_file" ]]; then
        echo "WARNING: Cannot write to log file: ${log_file}"
    fi
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

#' Convenience Logging Functions
for level in ${PROJECT_CONFIG[LOG_LEVELS]}; do
    level_lower=$(echo "$level" | tr '[:upper:]' '[:lower:]')
    eval "log_${level_lower}() { log_message \"$level\" \"\$1\" \"\$2\"; }"
done
