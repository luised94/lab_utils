#!/bin/bash
source "$HOME/lab_utils/bash/functions/filesystem_utils.sh"

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

validate_bmc_paths() {
    local bmc_server="$1"
    local experiment_id="$2"
    local log_file="$3"

    local bmc_path="${BMC_CONFIG[SOURCE_FS]}/$experiment_id"
    local local_path="$HOME/data/$experiment_id/fastq"
    
    # Add filesystem checks
    verify_filesystem_path "$bmc_path" "/net/bmc-pub17" "$log_file" || return 1
    check_filesystem_space "$local_path" "$log_file" || return 1
    
    mkdir -p "$local_path"
    echo -n "$bmc_path:$local_path"
}

download_from_bmc() {
    local paths="$1"
    local log_file="$2"

    local bmc_path=${paths%:*}
    local local_path=${paths#*:}
    
    log_info "Starting download from: $bmc_path" "$log_file"
    
    if ! srun rsync ${PROJECT_CONFIG[RSYNC_OPTIONS]} "$bmc_path/" "$local_path/"; then
        log_error "Download failed" "$log_file"
        return 1
    fi
    return 0
}
