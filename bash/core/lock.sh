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
