# bash/functions/lock_utils.sh
#!/bin/bash

source "$HOME/lab_utils/bash/config/logging_config.sh"

# Shared validation for both acquire and release

validate_lock_path() {
    local path="$1"
    local normalized_path
    
    # Basic checks
    [[ -z "$path" ]] && { echo "Empty path" >&2; return 1; }
    [[ "$path" == "/" ]] && { echo "Root path not allowed" >&2; return 1; }
    
    # Resolve path
    normalized_path=$(readlink -f "$path" 2>/dev/null) || {
        echo "Cannot resolve path: $path" >&2
        return 1
    }
    
    # Check against protected patterns
    for pattern in "${!PROTECTED_PATHS[@]}"; do
        if [[ "$normalized_path" =~ $pattern ]]; then
            echo "Protected path matched pattern: $pattern" >&2
            return 1
        fi
    done
    
    # Verify parent directory is writable
    local parent_dir=$(dirname "$normalized_path")
    if [[ ! -w "$parent_dir" ]]; then
        echo "Parent directory not writable: $parent_dir" >&2
        return 1
    fi
    
    return 0
}


#' Enhanced Lock Acquisition with Deadlock Prevention
acquire_lock() {
    local lock_file="$1"
    local timeout="${2:-${CORE_CONFIG[LOCK_TIMEOUT]}}"
    local start_time=$(date +%s)
    
    # Ensure clean state
    cleanup_stale_locks
    
    while true; do
        # Check timeout
        if (( $(date +%s) - start_time >= timeout )); then
            echo "ERROR: Lock acquisition timeout: $lock_file" >&2
            return 1
        }
        
        # Try to acquire lock
        if mkdir "$lock_file" 2>/dev/null; then
            # Record lock metadata
            echo $$ > "$lock_file/pid"
            date +%s > "$lock_file/timestamp"
            echo "${BASH_SOURCE[1]}" > "$lock_file/source"
            return 0
        fi
        
        # Check for stale lock
        if is_lock_stale "$lock_file"; then
            cleanup_lock "$lock_file"
            continue
        fi
        
        sleep 0.1
    done
}

#' Cleanup Stale Locks
cleanup_stale_locks() {
    local lock_dir="${CORE_CONFIG[LOCK_BASE_DIR]}"
    local max_age=3600  # 1 hour
    
    find "$lock_dir" -type d -name "*.lock" | while read -r lock; do
        is_lock_stale "$lock" && cleanup_lock "$lock"
    done
}

# Enhanced release_lock with shared validation
release_lock() {
    local lock_file="$1"
    local force="${2:-false}"
    
    if ! validate_lock_path "$lock_file"; then
        echo "ERROR: Invalid or protected lock path: $lock_file" >&2
        return 1
    fi
    
    [[ ! -d "$lock_file" ]] && return 0
    
    local pid_file="$lock_file/pid"
    if [[ -f "$pid_file" ]]; then
        local stored_pid=$(cat "$pid_file")
        if [[ "$stored_pid" != "$$" && "$force" != "true" ]]; then
            echo "ERROR: PID mismatch. Lock owned by $stored_pid" >&2
            return 1
        fi
    fi
    
    rm -rf "${lock_file:?}"/* 2>/dev/null
    rmdir "${lock_file:?}" 2>/dev/null
    return $?
}
