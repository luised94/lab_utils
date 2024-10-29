#!/bin/bash
# bash/scripts/download_bmc_data.sh

# Source dependencies
source "$HOME/lab_utils/bash/config/project_config.sh"
source "$HOME/lab_utils/bash/functions/logging_utils.sh"

#' Download BMC Data Main Function
#' @param bmc_server Character BMC server name
#' @param experiment_id Character Experiment identifier
#' @return Integer 0 if successful, 1 otherwise
download_bmc_data_main() {
    local log_file=$(initialize_logging "download_bmc_data")
    
    # Validate arguments
    if [[ $# -ne 2 ]]; then
        show_usage
        return 1
    fi
    local bmc_server="$1"
    local experiment_id="$2"
    
    # Validate paths
    local paths
    if ! paths=$(validate_bmc_paths "$bmc_server" "$experiment_id" "$log_file"); then
        return 1
    fi
    
    # Download data
    if ! download_from_bmc "$paths" "$log_file"; then
        return 1
    fi
    
    # Organize files
    if ! organize_fastq_files "${paths#*:}" "$log_file"; then
        return 1
    fi
    
    # Cleanup
    if ! cleanup_downloaded_data "${paths#*:}" "$log_file"; then
        return 1
    fi
    
    log_info "Download process completed successfully" "$log_file"
    return 0
}

#' Validate BMC Paths
#' @param bmc_server Character BMC server name
#' @param experiment_id Character Experiment identifier
#' @param log_file Character Log file path
#' @return String Combined paths or 1 if validation fails
validate_bmc_paths() {
    local bmc_server="$1"
    local experiment_id="$2"
    local log_file="$3"
    
    # Format BMC path
    local bmc_path=$(printf "${PROJECT_CONFIG[BMC_BASE_PATH]}" "$bmc_server" "$experiment_id")
    local local_path="${PROJECT_CONFIG[REMOTE_PATH]}/$experiment_id/${PROJECT_CONFIG[BMC_FASTQ_DIR]}"
    
    # Validate paths
    if [[ ! -d "$bmc_path" ]]; then
        log_error "BMC directory not found: $bmc_path" "$log_file"
        return 1
    fi
    
    if [[ ! -d "$local_path" ]]; then
        log_error "Local directory not found: $local_path" "$log_file"
        return 1
    fi
    
    echo "$bmc_path:$local_path"
}

#' Download Data from BMC
#' @param paths String Combined source:destination paths
#' @param log_file Character Log file path
#' @return Integer 0 if successful, 1 otherwise
download_from_bmc() {
    local paths="$1"
    local log_file="$2"
    
    local bmc_path=${paths%:*}
    local local_path=${paths#*:}
    
    log_info "Starting download from: $bmc_path" "$log_file"
    
    if ! srun rsync ${PROJECT_CONFIG[RSYNC_OPTIONS]} \
                    ${PROJECT_CONFIG[RSYNC_INCLUDES]} \
                    ${PROJECT_CONFIG[RSYNC_EXCLUDES]} \
                    "$bmc_path/" "$local_path/"; then
        log_error "Download failed" "$log_file"
        return 1
    fi
    
    return 0
}

# Show usage information
show_usage() {
    cat << EOF
Usage: $(basename "$0") <bmc_server> <experiment_id>

Arguments:
    bmc_server    BMC server name (e.g., ${PROJECT_CONFIG[BMC_DEFAULT_SERVER]})
    experiment_id Experiment identifier (format: YYMMDD'Bel')

Example:
    $(basename "$0") ${PROJECT_CONFIG[BMC_DEFAULT_SERVER]} 240808Bel
EOF
}

# Execute if run as script
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    download_bmc_data_main "$@"
fi
