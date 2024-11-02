#!/bin/bash
# bash/modules/fastq/bmc_handler.sh

#' BMC Data Handler Functions
#' @description Functions for BMC data management

source "${LAB_UTILS_ROOT}/bash/config/modules/bmc_config.sh"

verify_host() {
    local current_host=$(hostname)
    [[ "$current_host" == "luria" ]] || {
        log_error "Must run on luria.mit.edu (current: $current_host)"
        return 1
    }
}

validate_bmc_paths() {
    local experiment_id="$1"
    local log_file="$2"

    # Construct paths using exact locations
    local bmc_path="${BMC_CONFIG[SOURCE_FS]}/$experiment_id"
    local local_path="${BMC_CONFIG[TARGET_FS]}/$experiment_id/fastq"
    
    # Check source exists
    if [[ ! -d "$bmc_path" ]]; then
        log_error "BMC directory not found: $bmc_path" "$log_file"
        log_error "Please verify the experiment ID and try again" "$log_file"
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
