#!/bin/bash
#' Move fastq files inside directories to current working directory.
#' @param target_dir Character Target directory
#' @param log_file Character Log file path
#' @return Integer 0 if successful, 1 otherwise
move_fastq_files_to_current_directory() {
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
