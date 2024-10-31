#!/bin/bash
source "$HOME/lab_utils/bash/functions/filesystem_utils.sh"

remove_files_safely() {
    local type="$1"
    local pattern="$2"
    local log_file="$3"
    
    # Verify filesystem first
    local real_path=$(readlink -f .)
    verify_filesystem_path "$real_path" "${BMC_CONFIG[TARGET_FS]}" "$log_file" || return 1
    
    log_info "Searching for ${type}s matching: $pattern" "$log_file"
    
    # Find and sort items by size
    local items
    if [[ "$type" == "dir" ]]; then
        items=$(find . -type d -name "$pattern" -exec du -s {} \; | \
                sort -rn | \
                awk '{printf "%s\t%.2fGB\n", $2, $1/1024/1024}')
    else
        items=$(find . -type f -name "$pattern" -exec du -s {} \; | \
                sort -rn | \
                awk '{printf "%s\t%.2fGB\n", $2, $1/1024/1024}')
    fi
    
    # Check if items found
    if [[ -z "$items" ]]; then
        log_info "No ${type}s found matching: $pattern" "$log_file"
        return 0
    fi
    
    # Log items to be removed
    log_warning "Found ${type}s to remove:" "$log_file"
    echo "$items" | tee -a "$log_file"
    
    # Interactive confirmation if running interactively
    if [[ -t 0 ]]; then
        read -p "DELETE these ${type}s? (y/n) " -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            log_info "Deletion aborted by user" "$log_file"
            return 0
        fi
    fi
    
    # Remove items with size logging
    echo "$items" | while IFS=$'\t' read -r item size; do
        if [[ -n "$item" ]]; then
            log_info "Removing: $item (Size: $size)" "$log_file"
            if rm -rf "$item"; then
                log_info "Successfully removed: $item" "$log_file"
            else
                log_error "Failed to remove: $item" "$log_file"
            fi
        fi
    done
    
    return 0
}

cleanup_downloaded_data() {
    local target_dir="$1"
    local log_file="$2"
    
    verify_filesystem_path "$target_dir" "/net/bmc-pub14/data" "$log_file" || return 1
    
    log_info "Starting cleanup process in: $target_dir" "$log_file"
    
    # Store current directory
    local current_dir=$(pwd)
    
    # Change to target directory
    cd "$target_dir" || {
        log_error "Failed to access directory: $target_dir" "$log_file"
        return 1
    }
    
    # Process directory patterns
    IFS=' ' read -r -a dir_patterns <<< "${PROJECT_CONFIG[CLEANUP_DIRS]}"
    for pattern in "${dir_patterns[@]}"; do
        remove_files_safely "dir" "$pattern" "$log_file"
    done
    
    # Process file patterns
    IFS=' ' read -r -a file_patterns <<< "${PROJECT_CONFIG[CLEANUP_FILES]}"
    for pattern in "${file_patterns[@]}"; do
        remove_files_safely "file" "$pattern" "$log_file"
    done
    
    # Return to original directory
    cd "$current_dir"
    
    log_info "Cleanup process completed" "$log_file"
    return 0
}
