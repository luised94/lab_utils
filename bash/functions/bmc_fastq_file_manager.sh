#!/bin/bash
# bash/functions/fastq_file_manager.sh

#' Organize FASTQ Files
#' @param target_dir Character Target directory
#' @param log_file Character Log file path
#' @return Integer 0 if successful, 1 otherwise
organize_fastq_files() {
    local target_dir="$1"
    local log_file="$2"
    
    log_info "Organizing FASTQ files in: $target_dir" "$log_file"
    
    # Store current directory
    local current_dir=$(pwd)
    
    # Change to target directory
    cd "$target_dir" || {
        log_error "Failed to access directory: $target_dir" "$log_file"
        return 1
    }
    
    # Move files to root of target directory
    find . -type f -name "${PROJECT_CONFIG[FASTQ_PATTERN]}" -exec mv {} . \; || {
        log_error "Failed to move FASTQ files" "$log_file"
        cd "$current_dir"
        return 1
    }
    
    # Return to original directory
    cd "$current_dir"
    
    log_info "FASTQ files organized successfully" "$log_file"
    return 0
}

#' Safe Remove Function
#' @param type Character 'dir' or 'file'
#' @param pattern Character Pattern to match
#' @param log_file Character Log file path
#' @return Integer 0 if successful, 1 otherwise
safe_remove() {
    local type="$1"
    local pattern="$2"
    local log_file="$3"
    local staging_dir="/tmp/bmc_staging_$$"
    
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
    }
    
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

#' Cleanup Downloaded Data
#' @param target_dir Character Target directory
#' @param log_file Character Log file path
#' @return Integer 0 if successful, 1 otherwise
cleanup_downloaded_data() {
    local target_dir="$1"
    local log_file="$2"
    
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
