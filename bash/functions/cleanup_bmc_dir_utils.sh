#!/bin/bash
source "$HOME/lab_utils/bash/functions/filesystem_utils.sh"

safe_remove() {
    # Keep original implementation but add filesystem checks
    local type="$1"
    local pattern="$2"
    local log_file="$3"
    
    local real_path=$(readlink -f .)
    verify_filesystem_path "$real_path" "/net/bmc-pub14/data" "$log_file" || return 1
    
    log_info "Searching for ${type}s matching: $pattern" "$log_file"
    
    # Find matching items
    local items
    if [[ "$type" == "dir" ]]; then
        items=$(find . -type d -name "$pattern")
    else
        items=$(ls $pattern 2>/dev/null)
    fi
    
    # Check if items found
    if [[ -z "$items" ]]; then
        log_info "No ${type}s found matching: $pattern" "$log_file"
        return 0
    fi
    
    # Log items to be removed
    log_warning "Found ${type}s to remove:" "$log_file"
    echo "$items" | tee -a "$log_file"
    
    # Interactive confirmation
    if [[ -t 0 ]]; then  # Check if running interactively
        read -p "DELETE these ${type}s? (y/n) " -r
        echo
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            log_info "Deletion aborted by user" "$log_file"
            return 0
        fi
    fi
    
    # Create staging directory
    mkdir -p "$staging_dir" || {
        log_error "Failed to create staging directory" "$log_file"
        return 1
    }
    
    # Move and remove items
    echo "$items" | while read item; do
        if [[ -n "$item" ]]; then
            mv "$item" "$staging_dir/" || {
                log_warning "Failed to move: $item" "$log_file"
            }
        fi
    done
    
    # Remove staging directory
    rm -rf "$staging_dir"
    log_info "Cleanup completed successfully" "$log_file"
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
        safe_remove "dir" "$pattern" "$log_file"
    done
    
    # Process file patterns
    IFS=' ' read -r -a file_patterns <<< "${PROJECT_CONFIG[CLEANUP_FILES]}"
    for pattern in "${file_patterns[@]}"; do
        safe_remove "file" "$pattern" "$log_file"
    done
    
    # Return to original directory
    cd "$current_dir"
    
    log_info "Cleanup process completed" "$log_file"
    return 0
}
