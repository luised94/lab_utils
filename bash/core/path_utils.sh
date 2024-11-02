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
