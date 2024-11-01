# bash/functions/lock_utils.sh
#!/bin/bash

source "$HOME/lab_utils/bash/config/lock_config.sh"

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
    }
    
    return 0
}

# Modified acquire_lock to use shared validation
acquire_lock() {
    local lock_file="$1"
    local timeout="${2:-$LOGGING_CONFIG[LOCK_TIMEOUT]}"
    
    if ! validate_lock_path "$lock_file"; then
        echo "ERROR: Invalid or protected lock path: $lock_file" >&2
        return 1
    }
    
    local wait_time=0
    while [[ $wait_time -lt $timeout ]]; do
        if mkdir "$lock_file" 2>/dev/null; then
            echo $$ > "$lock_file/pid"
            return 0
        fi
        sleep 1
        ((wait_time++))
    done
    return 1
}

# Enhanced release_lock with shared validation
release_lock() {
    local lock_file="$1"
    local force="${2:-false}"
    
    if ! validate_lock_path "$lock_file"; then
        echo "ERROR: Invalid or protected lock path: $lock_file" >&2
        return 1
    }
    
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
