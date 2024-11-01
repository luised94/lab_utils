# bash/functions/lock_utils.sh
#!/bin/bash

#' Acquire File Lock
#' @param lock_file Character Path to lock file
#' @param timeout Integer Seconds to wait
#' @return Integer 0 if successful
acquire_lock() {
    local lock_file="$1"
    local timeout="${2:-${LOGGING_CONFIG[LOCK_TIMEOUT]}}"
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

#' Release File Lock
#' @param lock_file Character Path to lock file
release_lock() {
    local lock_file="$1"
    rm -rf "$lock_file" 2>/dev/null
}
