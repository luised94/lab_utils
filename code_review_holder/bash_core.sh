
#!/bin/bash
source "$HOME/lab_utils/bash/config/bmc_config.sh"

remove_files_safely() {
    local pattern="$1"
    local log_file="$2"
    
    log_info "Removing items matching: $pattern" "$log_file"
    
    # Find matching items with sizes
    local items=$(find . -name "$pattern" -exec du -sh {} \; | \
                 sort -hr)
    
    if [[ -z "$items" ]]; then
        log_info "No items found matching: $pattern" "$log_file"
        return 0
    fi
    
    # Log items to be removed
    log_warning "Will remove:" "$log_file"
    echo "$items" | tee -a "$log_file"
    
    # Interactive confirmation if running interactively
    if [[ -t 0 ]]; then
        read -p "Proceed with deletion? (y/n) " -r
        [[ ! $REPLY =~ ^[Yy]$ ]] && return 0
    fi
    
    # Remove items
    echo "$items" | while read size item; do
        if rm -rf "${item#./}"; then
            log_info "Removed: $item ($size)" "$log_file"
        else
            log_error "Failed to remove: $item" "$log_file"
        fi
    done
}

#' Check Available Filesystem Space
#' @param path Character Path to check
#' @param log_file Character Log file path
#' @param min_space Numeric Minimum required space in GB
#' @return Integer 0 if sufficient space
check_filesystem_space() {
    local path="$1"
    local log_file="$2"
    local min_space="${3:-${CORE_CONFIG[MIN_SPACE_GB]:-10}}"  # Default 10GB if not set
    
    # Move to CORE_CONFIG
    # [MIN_SPACE_GB]="10"
    
    local available_gb=$(df -P "$path" | awk 'NR==2 {print $4/1024/1024}')
    
    log_info "Available space: ${available_gb}GB" "$log_file"
    
    if (( $(echo "$available_gb < $min_space" | bc -l) )); then
        log_error "Insufficient space: need ${min_space}GB, have ${available_gb}GB" "$log_file"
        return 1
    fi
    return 0
}

#' Verify Path is Under Expected Filesystem
#' @param path Character Path to verify
#' @param expected_fs Character Expected filesystem root
#' @param log_file Character Log file path
#' @return Integer 0 if path is valid
verify_filesystem_path() {
    local path="$1"
    local expected_fs="$2"
    local log_file="$3"
    
    local real_path=$(readlink -f "$path")
    if [[ "$real_path" != "$expected_fs"* ]]; then
        log_error "Path not under expected filesystem: $path" "$log_file"
        log_error "Expected root: $expected_fs" "$log_file"
        return 1
    fi
    return 0
}


#!/bin/bash

# Load required settings for logging_utils
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
    local retry_count=0

    # Ensure lock directory exists
    mkdir -p "${CORE_CONFIG[LOCK_BASE_DIR]}" 2>/dev/null
    mkdir -p "${CORE_CONFIG[DEFAULT_LOG_ROOT]}" 2>/dev/null

    while [[ $retry_count -lt ${CORE_CONFIG[LOCK_RETRY]} ]]; do
        if acquire_lock "$lock_name"; then
            echo "$entry" >> "$log_file"
            local status=$?
            release_lock "$lock_name"
            return $status
        fi
        ((retry_count++))
        sleep 1
    done
    
    echo "ERROR: Failed to acquire log lock after ${CORE_CONFIG[LOCK_RETRY]} attempts" >&2
    return 1
}

#' Get Run Count
#' @param log_file Character Path to log file
#' @return Integer Run count
get_run_count() {
    local log_file="$1"
    local lock_name="log_$(basename "$log_file")"
    local count=0
    
    # Acquire lock for counting
    if ! acquire_lock "$lock_name"; then
        echo "0"
        return 1
    fi
    
    # Count runs with error handling
    if [[ -f "$log_file" ]]; then
        count=$(grep -c "${CORE_CONFIG[RUN_SEPARATOR]}" "$log_file" 2>/dev/null || echo "0")
    fi
    
    release_lock "$lock_name"
    echo "$count"
}

#' Format Run Entry
#' @param count Integer Run count
#' @return String Formatted entry
format_run_entry() {
    local count="$1"
    local timestamp
    timestamp=$(date +"${CORE_CONFIG[TIMESTAMP_FORMAT]}")
    
    if [[ $count -eq 0 ]]; then
        printf "${CORE_CONFIG[FIRST_RUN_FORMAT]}" \
            "${CORE_CONFIG[RUN_SEPARATOR]}" \
            "$timestamp"
    else
        printf "${CORE_CONFIG[ENTRY_FORMAT]}" \
            "${CORE_CONFIG[RUN_SEPARATOR]}" \
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
    local log_dir="${2:-${CORE_CONFIG[DEFAULT_LOG_ROOT]}}"
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
    echo "$log_file"
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
    
    # Improved level validation
    local valid_levels=(${CORE_CONFIG[LOG_LEVELS]})
    local is_valid=0
    for valid_level in "${valid_levels[@]}"; do
        if [[ "$level" == "$valid_level" ]]; then
            is_valid=1
            break
        fi
    done
    
    if [[ $is_valid -eq 0 ]]; then
        echo "Invalid log level: $level (Valid levels: ${CORE_CONFIG[LOG_LEVELS]})" >&2
        return 1
    fi
    
    # Truncate long messages
    if [[ ${#message} -gt ${CORE_CONFIG[MAX_MESSAGE_LENGTH]} ]]; then
        message="${message:0:${CORE_CONFIG[MAX_MESSAGE_LENGTH]}}..."
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
for level in ${CORE_CONFIG[LOG_LEVELS]}; do
    level_lower=$(echo "$level" | tr '[:upper:]' '[:lower:]')
    eval "log_${level_lower}() { log_message \"$level\" \"\$1\" \"\$2\"; }"
done

#!/bin/bash
# bash/core/initialize_lab_environment.sh

#' Initialize Lab Utils Environment
#' @description Core initialization script for lab utilities
#' @export LAB_UTILS_ROOT, LAB_UTILS_INITIALIZED

# Guard against multiple inclusion
[[ -n "$LAB_UTILS_INITIALIZED" ]] && return

#' Discover Lab Utils Root Using Git
#' @description Find repository root using git
#' @return String Absolute path to repository root
discover_lab_utils_root() {
    local repo_root
    
    # Check if we're in a git repository
    if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
        echo "? Not in a git repository" >&2
        return 1
    fi
    
    # Get repository root (compatible with git 1.8+)
    repo_root="$(git rev-parse --show-toplevel 2>/dev/null)" || {
        echo "? Failed to determine repository root" >&2
        return 1
    }
    
    # Verify it's the correct repository
    if [[ ! -d "$repo_root/bash/core" ]]; then
        echo "? Invalid repository structure" >&2
        return 1
    fi
    
    echo "$repo_root"
}

# bash/core/initialize_lab_environment.sh

#' Load Module Configuration
#' @param module_name Character Name of module
#' @return Integer 0 if successful
load_module_config() {
    local module_name="$1"
    local config_path="${LAB_UTILS_ROOT}/bash/config/modules/${module_name}_config.sh"
    
    log_debug "Loading module configuration: ${module_name}" "${LAB_UTILS_LOG_FILE:-/dev/null}"
    
    if [[ ! -f "$config_path" ]]; then
        log_error "Module configuration not found: ${module_name}_config.sh"
        return 1
    fi
    
    source "$config_path" || {
        log_error "Failed to load module configuration: ${module_name}"
        return 1
    }
    
    return 0
}

#' Load Laboratory Module
#' @param module_path Character Path relative to modules directory
#' @return Integer 0 if successful
load_lab_module() {
    local module_path="$1"
    local module_name="${module_path%%/*}"  # Extract base module name
    
    # Load module configuration first
    load_module_config "$module_name" || {
        log_error "Failed to load configuration for module: $module_name"
        return 1
    }
    
    # Then load module implementation
    local full_path="${LAB_UTILS_ROOT}/bash/modules/${module_path}.sh"
    if [[ ! -f "$full_path" ]]; then
        log_error "Module not found: $module_path"
        return 1
    fi
    
    source "$full_path" || {
        log_error "Failed to load module: $module_path"
        return 1
    }
    
    return 0
}
#' Load Laboratory Module
#' @description Load module from lab utils repository
#' @param module_path Character Path relative to modules directory
#' @return Integer 0 if successful
load_lab_module() {
    local module_path="$1"
    local full_path="${LAB_UTILS_ROOT}/bash/modules/${module_path}.sh"
    
    # Debug output
    log_debug "Loading module: $module_path" "${LAB_UTILS_LOG_FILE:-/dev/null}"
    
    # Validate path
    if [[ ! -f "$full_path" ]]; then
        log_error "Module not found: $module_path" "${LAB_UTILS_LOG_FILE:-/dev/null}"
        return 1
    fi
    
    # Source module with error handling
    if ! source "$full_path"; then
        log_error "Failed to load module: $module_path" "${LAB_UTILS_LOG_FILE:-/dev/null}"
        return 1
    fi
    
    log_debug "Successfully loaded: $module_path" "${LAB_UTILS_LOG_FILE:-/dev/null}"
    return 0
}

# Set root directory
if [[ -z "$LAB_UTILS_ROOT" ]]; then
    LAB_UTILS_ROOT="$(discover_lab_utils_root)"
    echo "DEBUG: LAB_UTILS_ROOT set to: $LAB_UTILS_ROOT"
    echo "DEBUG: Looking for config in: $LAB_UTILS_ROOT/bash/config"
    export LAB_UTILS_ROOT
fi

# Verify critical directory exists
if [[ ! -d "$LAB_UTILS_ROOT" ]]; then
    echo "ERROR: Lab utils directory not found: $LAB_UTILS_ROOT" >&2
    return 1
fi

# Core configuration files (order matters)
readonly CORE_CONFIG_FILES=(
    "core_config.sh"      # Base configuration with logging and lock settings.
)

# Core module files (order matters)
readonly CORE_MODULES=(
    "logging.sh"          # Must be first
    "lock.sh"
    "path_utils.sh"
)

# Initialize core configuration
for config in "${CORE_CONFIG_FILES[@]}"; do
    config_path="${LAB_UTILS_ROOT}/bash/config/${config}"
    if [[ -f "$config_path" ]]; then
        source "$config_path"
    else
        echo "ERROR: Required configuration not found: ${config}" >&2
        return 1
    fi
done

# Initialize core modules
for module in "${CORE_MODULES[@]}"; do
    module_path="${LAB_UTILS_ROOT}/bash/core/${module}"
    if [[ -f "$module_path" ]]; then
        source "$module_path"
    else
        echo "ERROR: Required module not found: ${module}" >&2
        return 1
    fi
done

# Export common paths
export LAB_UTILS_CONFIG_DIR="${LAB_UTILS_ROOT}/bash/config"
export LAB_UTILS_CORE_DIR="${LAB_UTILS_ROOT}/bash/core"
export LAB_UTILS_MODULES_DIR="${LAB_UTILS_ROOT}/bash/modules"
export LAB_UTILS_SCRIPTS_DIR="${LAB_UTILS_ROOT}/bash/scripts"
export LAB_UTILS_TESTS_DIR="${LAB_UTILS_ROOT}/bash/tests"


source "${LAB_UTILS_CORE_DIR}/config_export.sh" || {
    echo "[ERROR] Failed to load config export functions" >&2
    return 1
}

# After loading config
export_core_config || {
    echo "[ERROR] Failed to export configuration" >&2
    return 1
}

# Verify export worked
if [[ -z "$LAB_UTILS_CONFIG_SERIALIZED" ]]; then
    echo "[ERROR] Configuration export failed" >&2
    return 1
fi
# Export critical functions
export -f import_core_config
export -f export_core_config

# Initialize logging
if [[ -z "$LAB_UTILS_LOG_FILE" ]]; then
    initialize_logging "lab_utils_init"
fi

# Cleanup stale locks on initialization
#cleanup_stale_locks 2>/dev/null

# Mark as initialized
readonly LAB_UTILS_INITIALIZED=1

log_info "Lab environment initialized successfully"

#!/bin/bash
# bash/core/config_export.sh

export_core_config() {
    local serialized=""
    for key in "${!CORE_CONFIG[@]}"; do
        serialized+="${key}=${CORE_CONFIG[$key]};"
    done
    export LAB_UTILS_CONFIG_SERIALIZED="${serialized}"
    
    # Export individual keys
    for key in "${!CORE_CONFIG[@]}"; do
        export "LAB_UTILS_${key}=${CORE_CONFIG[$key]}"
    done
}

import_core_config() {
    declare -g -A CORE_CONFIG
    while IFS='=' read -r key value; do
        [[ -n "$key" ]] && CORE_CONFIG[$key]="$value"
    done < <(echo "$LAB_UTILS_CONFIG_SERIALIZED" | tr ';' '\n')
}

#!/bin/bash
# bash/functions/lock_utils.sh

# Source core_config and logging or through initialization.
#source "$HOME/lab_utils/bash/config/logging_config.sh"

# Shared validation for both acquire and release
#' Validate Lock Path
#' @param lock_path Character Path to lock file
#' @return Integer 0 if valid
validate_lock_path() {
    local lock_path="$1"
    local normalized_path
    
    # Basic checks
    [[ -z "$lock_path" ]] && { 
        log_error "Empty lock path"
        return 1
    }
    
    # Ensure path is under user's lock directory
    normalized_path=$(readlink -f "$lock_path" 2>/dev/null) || {
        log_error "Cannot resolve path: $lock_path"
        return 1
    }
    
    local user_lock_dir="${CORE_CONFIG[LOCK_BASE_DIR]}/${USER}"
    if [[ ! "$normalized_path" =~ ^"$user_lock_dir" ]]; then
        log_error "Lock must be in user directory: $user_lock_dir"
        return 1
    fi
    
    # Check parent directory permissions
    local parent_dir=$(dirname "$normalized_path")
    if [[ ! -w "$parent_dir" ]]; then
        log_error "Parent directory not writable: $parent_dir"
        return 1
    fi
    
    return 0
}

#' Acquire Lock
#' @param lock_name Character Lock identifier
#' @param timeout Integer Seconds to wait
acquire_lock() {
    local lock_name="$1"
    local timeout="${2:-${CORE_CONFIG[LOCK_TIMEOUT]}}"
    local user_specific_dir="${CORE_CONFIG[LOCK_BASE_DIR]}/${USER}"
    local lock_file="$user_specific_dir/${lock_name}.lock"
    
    # Ensure user lock directory exists
    create_lock_directory
    
    # Attempt to acquire lock
    local start_time=$(date +%s)
    while (( $(date +%s) - start_time < timeout )); do
        if mkdir "$lock_file" 2>/dev/null; then
            echo $$ > "$lock_file/pid"
            chmod 700 "$lock_file"  # Only owner can access
            return 0
        fi
        sleep 0.1
    done
    
    return 1
}

#' Cleanup Stale Locks
cleanup_stale_locks() {
    local lock_dir="${CORE_CONFIG[LOCK_BASE_DIR]}"
    local max_age=3600  # 1 hour
    local user_specific_dir="${lock_dir}/${USER}"
    
    # Only clean user-specific locks
    if [[ ! -d "$user_specific_dir" ]]; then
        return 0
    fi
    
    # Find only locks owned by current user
    find "$user_specific_dir" -type d -name "*.lock" -user "$USER" -mmin +60 | \
    while read -r lock; do
        if validate_lock_path "$lock"; then
            local pid_file="$lock/pid"
            if [[ -f "$pid_file" ]]; then
                local stored_pid=$(cat "$pid_file")
                # Check if process still exists
                if ! kill -0 "$stored_pid" 2>/dev/null; then
                    rm -rf "$lock"
                    log_debug "Removed stale lock: $lock (PID: $stored_pid)"
                fi
            else
                rm -rf "$lock"
                log_debug "Removed invalid lock: $lock (no PID file)"
            fi
        fi
    done
    
    log_info "Cleaned up stale locks for user: $USER"
}

#' Create Lock Directory
#' @description Ensure user-specific lock directory exists
create_lock_directory() {
    local lock_base="${CORE_CONFIG[LOCK_BASE_DIR]}"
    local user_dir="$lock_base/$USER"
    
    # Create with strict permissions
    mkdir -p "$lock_base"
    chmod 1777 "$lock_base"  # Like /tmp
    
    mkdir -p "$user_dir"
    chmod 700 "$user_dir"    # Only user can access
    
    return 0
}

#' Release Lock
#' @param lock_name Character Lock identifier
#' @param force Logical Force release even if not owner
#' @return Integer 0 if successful
release_lock() {
    local lock_name="$1"
    local force="${2:-false}"
    local user_specific_dir="${CORE_CONFIG[LOCK_BASE_DIR]}/${USER}"
    local lock_file="$user_specific_dir/${lock_name}.lock"
    local verbose="${3:-${CORE_CONFIG[VERBOSE]:-false}}"  # ÃÄ Add verbose parameter
    
    # Validate lock path
    if ! validate_lock_path "$lock_file"; then
        log_error "Invalid or protected lock path: $lock_file"
        return 1
    fi
    
    # Check if lock exists
    if [[ ! -d "$lock_file" ]]; then
        log_debug "Lock already released: $lock_file"
        return 0
    fi
    
    # Verify ownership
    if [[ ! -O "$lock_file" ]]; then
        log_error "Lock owned by different user: $(stat -c %U "$lock_file")"
        return 1
    fi
    
    # Check PID
    local pid_file="$lock_file/pid"
    if [[ -f "$pid_file" ]]; then
        local stored_pid
        stored_pid=$(cat "$pid_file" 2>/dev/null)
        
        if [[ "$stored_pid" != "$$" && "$force" != "true" ]]; then
            # Additional check: see if process is still running
            if kill -0 "$stored_pid" 2>/dev/null; then
                log_error "Lock owned by running process: $stored_pid"
                return 1
            elif [[ "$force" != "true" ]]; then
                log_warning "Found stale lock from dead process: $stored_pid"
            fi
        fi
    fi
    
    # Safe removal with error checking
    {
        rm -f "${lock_file:?}/pid" 2>/dev/null
        rm -rf "${lock_file:?}"/* 2>/dev/null
        rmdir "${lock_file:?}" 2>/dev/null
    } || {
        log_error "Failed to remove lock: $lock_file"
        return 1
    }

    if [[ "$verbose" == "true" ]]; then
        log_debug "Released lock: $lock_name"  # ÀÄ Only log if verbose
    fi
    return 0
}
