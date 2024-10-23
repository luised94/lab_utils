#!/bin/bash

source "../config/data_transfer_config.sh"

function validate_bmc_paths() {
    local server="$1"
    local experiment="$2"
    
    log_info "Validating BMC paths"
    
    local bmc_path=$(printf "${BMC_CONFIG[BASE_PATH]}" "$server" "$experiment")
    local local_path="${BMC_CONFIG[LOCAL_BASE]}/$experiment/${BMC_CONFIG[FASTQ_DIR]}"
    
    if [ ! -d "$bmc_path" ]; then
        log_error "BMC directory not found: $bmc_path"
        log_error "Verify BMC server and experiment name"
        return 1
    }
    
    if [ ! -d "$local_path" ]; then
        log_error "Local directory not found: $local_path"
        log_error "Run setup_experiment_dir first"
        return 1
    }
    
    echo "$bmc_path:$local_path"
}

function download_bmc_data() {
    local server="$1"
    local experiment="$2"
    
    local paths=$(validate_bmc_paths "$server" "$experiment") || return 1
    local bmc_path=${paths%:*}
    local local_path=${paths#*:}
    
    log_info "Downloading from BMC server: $server"
    
    local rsync_cmd="srun rsync ${RSYNC_CONFIG[OPTIONS]} \
                     ${RSYNC_CONFIG[INCLUDES]} \
                     ${RSYNC_CONFIG[EXCLUDES]} \
                     $bmc_path $local_path"
    
    if ! eval "$rsync_cmd"; then
        log_error "Rsync failed. Check connection or permissions"
        return 1
    }
    
    log_info "Download completed successfully"
    return 0
}
