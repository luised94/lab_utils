#!/bin/bash
source "$HOME/lab_utils/bash/functions/filesystem_utils.sh"

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
