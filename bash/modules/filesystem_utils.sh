#!/bin/bash
source "$HOME/lab_utils/bash/config/bmc_config.sh"
source "$HOME/lab_utils/bash/functions/logging_utils.sh"

check_filesystem_space() {
    local path="$1"
    local log_file="$2"
    local min_space=${3:-${BMC_CONFIG[MIN_SPACE_GB]}}
    
    local available_gb=$(df -P "$path" | awk 'NR==2 {print $4/1024/1024}')
    
    log_info "  Available space: ${available_gb}GB" "$log_file"
    
    if (( $(echo "$available_gb < $min_space" | bc -l) )); then
        log_error "Insufficient space, need ${min_space}GB, have ${available_gb}GB" "$log_file"
        return 1
    fi
    return 0
}

verify_filesystem_path() {
    local path="$1"
    local expected_fs="$2"
    local log_file="$3"
    
    local real_path=$(readlink -f "$path")
    if [[ "$real_path" != "$expected_fs"* ]]; then
        log_error "Invalid filesystem: $path" "$log_file"
        return 1
    fi
    return 0
}
