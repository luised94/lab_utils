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
    find . -type f -name "*.fastq" -exec mv {} . \; || {
        log_error "Failed to move FASTQ files" "$log_file"
        cd "$current_dir"
        return 1
    }
    # Return to original directory
    cd "$current_dir"
    log_info "FASTQ files organized successfully" "$log_file"
    return 0
}

clean_experiment_directory() {
    local target_dir="$1"
    local log_file="$2"
    log_info "Starting cleanup process in: $target_dir" "$log_file"
    # Store current directory
    local current_dir=$(pwd)
    # Verify we're in the correct directory structure
    if [[ "$target_dir" != "${BMC_CONFIG[TARGET_FS]}"* ]]; then
        log_error "Invalid target directory: $target_dir" "$log_file"
        log_error "Must be under: ${BMC_CONFIG[TARGET_FS]}" "$log_file"
        return 1
    fi
    # Change to target directory
    cd "$target_dir" || {
        log_error "Failed to access directory: $target_dir" "$log_file"
        return 1
    }

    # Process cleanup patterns
    local -a patterns=(${BMC_CONFIG[CLEANUP_DIRS]} ${BMC_CONFIG[CLEANUP_FILES]})
    for pattern in "${patterns[@]}"; do
        remove_files_safely "$pattern" "$log_file"
    done
    # Return to original directory
    cd "$current_dir"
    
    log_info "Cleanup process completed" "$log_file"
    return 0
}
