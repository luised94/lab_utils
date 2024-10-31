#!/bin/bash
source "$HOME/lab_utils/bash/functions/filesystem_utils.sh"

#' Verify Host Environment
#' @return Integer 0 if valid host, 1 otherwise
verify_host() {
    local current_host=$(hostname)
    if [[ "$current_host" != "luria" ]]; then
        log_error "This script must be run on luria.mit.edu"
        log_error "Current host: $current_host"
        return 1
    fi
    return 0
}


validate_bmc_paths() {
    local bmc_server="$1"
    local experiment_id="$2"
    local log_file="$3"

    # Construct paths using exact locations
    local bmc_path="${BMC_CONFIG[SOURCE_FS]}/$experiment_id"
    local local_path="${BMC_CONFIG[TARGET_FS]}/$experiment_id/fastq"
    
    log_debug "Validating paths:" "$log_file"
    log_debug "  BMC path: $bmc_path" "$log_file"
    log_debug "  Local path: $local_path" "$log_file"
    
    # Check source exists
    if [[ ! -d "$bmc_path" ]]; then
        log_error "BMC directory not found: $bmc_path" "$log_file"
        return 1
    fi
    
    # Check space on target filesystem
    if ! check_filesystem_space "$local_path" "$log_file"; then
        return 1
    fi
    
    # Create directory
    if ! mkdir -p "$local_path"; then
        log_error "Failed to create directory: $local_path" "$log_file"
        return 1
    fi
    
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
