
#!/bin/bash

# Load required settings for logging_utils
source "$HOME/lab_utils/bash/config/project_config.sh"
source "$HOME/lab_utils/bash/config/logging_config.sh"
source "$HOME/lab_utils/bash/functions/lock_utils.sh"
# Advanced Logging Functions for Bash Scripts
#
# Script: 002_logging_functions.sh
# Description: A set of functions for consistent logging across Bash scripts
# Author: Your Name
# Date: 2024-10-18

# Function: extract_script_path
# Purpose: Extract the full path of the current script
# Parameters: None
# Return: Script path
extract_script_path() {
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


#' Write Log Entry Atomically
#' @param entry Character Log entry
#' @param log_file Character Log file path
#' @return Integer 0 if successful
write_log_atomic() {
    local entry="$1"
    local log_file="$2"
    local lock_name="log_$(basename "$log_file")"
    local lock_file="${LOGGING_CONFIG[LOCK_BASE_DIR]}/${lock_name}.lock"
    local retry_count=0

    # Ensure lock directory exists
    mkdir -p "${LOGGING_CONFIG[LOCK_BASE_DIR]}" 2>/dev/null

    while [[ $retry_count -lt ${LOGGING_CONFIG[LOCK_RETRY]} ]]; do
        if acquire_lock "$lock_file"; then
            echo "$entry" >> "$log_file"
            local status=$?
            release_lock "$lock_file"
            return $status
        fi
        ((retry_count++))
        sleep 1
    done
    
    echo "ERROR: Failed to acquire log lock after ${LOGGING_CONFIG[LOCK_RETRY]} attempts" >&2
    return 1
}

#' Get Run Count
#' @param log_file Character Path to log file
#' @return Integer Run count
get_run_count() {
    local log_file="$1"
    local lock_file="${log_file}.lock"
    local count=0
    
    # Acquire lock for counting
    if ! acquire_lock "$lock_file"; then
        echo "0"
        return 1
    fi
    
    # Count runs with error handling
    if [[ -f "$log_file" ]]; then
        count=$(grep -c "${LOGGING_CONFIG[RUN_SEPARATOR]}" "$log_file" 2>/dev/null || echo "0")
    fi
    
    release_lock "$lock_file"
    echo "$count"
}

#' Format Run Entry
#' @param count Integer Run count
#' @return String Formatted entry
format_run_entry() {
    local count="$1"
    local timestamp
    timestamp=$(date +"${LOGGING_CONFIG[TIMESTAMP_FORMAT]}")
    
    if [[ $count -eq 0 ]]; then
        printf "${LOGGING_CONFIG[FIRST_RUN_FORMAT]}" \
            "${LOGGING_CONFIG[RUN_SEPARATOR]}" \
            "$timestamp"
    else
        printf "${LOGGING_CONFIG[ENTRY_FORMAT]}" \
            "${LOGGING_CONFIG[RUN_SEPARATOR]}" \
            "$((count + 1))" \
            "$timestamp"
    fi
}
#' Initialize Logging
#' @param script_name Character Name of the calling script
#' @param log_dir Character Optional log directory
#' @return String Path to log file

initialize_logging() {
    local script_name="${1:-$(basename "${BASH_SOURCE[1]}" )}"
    local log_dir="${2:-${LOGGING_CONFIG[DEFAULT_LOG_ROOT]}}"
    local log_file
    
    # Setup log file
    log_file="$log_dir/$(date +%Y-%m)/$(date +%Y-%m-%d)_${script_name}.log"
    mkdir -p "$(dirname "$log_file")"
    
    # Get run count atomically
    local run_count
    run_count=$(get_run_count "$log_file")
    
    # Write header atomically
    local entry
    entry=$(format_run_entry "$run_count")
    write_log_atomic "$entry" "$log_file"
    log_system_info "$log_file"
    log_git_info "$log_file"
    echo -n "$log_file"
}

#' Enhanced Log Message
#' @param level Character Log level
#' @param message Character Message to log
#' @param log_file Character Log file path
#' @return Integer 0 if successful
log_message() {
    local level="$1"
    local message="$2"
    local log_file="$3"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    local task_id="${SLURM_ARRAY_TASK_ID:-standalone}"
    local job_id="${SLURM_JOB_ID:-local}"
    # Validate level
    if [[ ! " ${LOGGING_CONFIG[LOG_LEVELS]} " =~ " ${level} " ]]; then
        echo "Invalid log level: $level" >&2
        return 1
    fi
    # Truncate long messages
    if [[ ${#message} -gt ${LOGGING_CONFIG[MAX_MESSAGE_LENGTH]} ]]; then
        message="${message:0:${LOGGING_CONFIG[MAX_MESSAGE_LENGTH]}}..."
    fi
    # Format entry
    local log_entry="[${timestamp}] [${level}] [Job:${job_id}] [Task:${task_id}] ${message}"
    # Console output
    echo "$log_entry" >&2
    # File output with atomic writes
    if [[ -n "$log_file" ]]; then
        write_log_atomic "$log_entry" "$log_file"
        return $?
    fi
    
    return 0
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
for level in ${LOGGING_CONFIG[LOG_LEVELS]}; do
    level_lower=$(echo "$level" | tr '[:upper:]' '[:lower:]')
    eval "log_${level_lower}() { log_message \"$level\" \"\$1\" \"\$2\"; }"
done
