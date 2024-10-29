#!/bin/bash
# bash/scripts/download_bmc_data.sh

# Source dependencies
source "$HOME/lab_utils/bash/config/project_config.sh"
source "$HOME/lab_utils/bash/functions/logging_utils.sh"
source "$HOME/lab_utils/bash/functions/bmc_fastq_file_manager.sh"

#' Verify Host Environment
#' @return Integer 0 if valid host, 1 otherwise
verify_host() {
    local current_host=$(hostname)
    if [[ "$current_host" != "luria" ]]; then
        echo "ERROR: This script must be run on luria.mit.edu"
        echo "Current host: $current_host"
        echo "Please transfer data to luria first and run this script there."
        return 1
    fi
    return 0
}

#' Download BMC Data Main Function
#' @param bmc_server Character BMC server name
#' @param experiment_id Character Experiment identifier
#' @return Integer 0 if successful, 1 otherwise
download_bmc_data_main() {

    # Verify host first
    if ! verify_host; then
        return 1
    fi

    # Validate arguments
    if [[ $# -ne 2 ]]; then
        show_usage
        log_error "Invalid number of arguments"
        return 1
    fi

    local bmc_server="$1"
    local experiment_id="$2"

        # Initialize logging with error checking
    local log_file
    log_file="$(initialize_logging "download_bmc_fastq_data")"
    if [[ -z "$log_file" ]]; then
        echo "ERROR: Failed to initialize logging"
        return 1
    fi
    log_info "Starting BMC data download for experiment: $experiment_id" "$log_file"
    echo "$log_file"
    echo "See line right above."

    # Verify log file exists and is writable
    if [[ ! -w "$log_file" ]]; then
        echo "ERROR: Log file not writable: $log_file"
        return 1
    fi


    log_trace "001: Made it here."
    echo "$log_file"
    # Validate paths
    local paths
    if ! paths=$(validate_bmc_paths "$bmc_server" "$experiment_id" "$log_file"); then
        return 1
    fi

    log_trace "002: After validate_bmc_paths"
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

#!/bin/bash
# Temporary debug version of validate_bmc_paths

validate_bmc_paths() {
    local bmc_server="$1"
    local experiment_id="$2"
    local log_file="$3"

    echo "=== Debug Information ==="
    echo "Arguments:"
    echo "  bmc_server: $bmc_server"
    echo "  experiment_id: $experiment_id"
    echo "  log_file: $log_file"

    echo "Configuration:"
    echo "  BMC_BASE_PATH: ${PROJECT_CONFIG[BMC_BASE_PATH]}"
    echo "  REMOTE_PATH: ${PROJECT_CONFIG[REMOTE_PATH]}"
    echo "  BMC_FASTQ_DIR: ${PROJECT_CONFIG[BMC_FASTQ_DIR]}"

    # Format BMC path
    local bmc_path
    printf -v bmc_path "${PROJECT_CONFIG[BMC_BASE_PATH]}" "$bmc_server" "$experiment_id"
    echo "Constructed paths:"
    echo "  BMC path: $bmc_path"

    local local_path="${PROJECT_CONFIG[REMOTE_PATH]}/$experiment_id/${PROJECT_CONFIG[BMC_FASTQ_DIR]}"
    echo "  Local path: $local_path"

    # Directory checks
    echo "Directory status:"
    echo "  BMC path exists: $([[ -d "$bmc_path" ]] && echo "yes" || echo "no")"
    echo "  Local path exists: $([[ -d "$local_path" ]] && echo "yes" || echo "no")"

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
