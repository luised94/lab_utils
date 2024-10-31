#!/bin/bash
# bash/scripts/download_bmc_data.sh

# Source dependencies
source "$HOME/lab_utils/bash/config/project_config.sh"
source "$HOME/lab_utils/bash/config/bmc_config.sh"
source "$HOME/lab_utils/bash/functions/filesystem_utils.sh"
source "$HOME/lab_utils/bash/functions/logging_utils.sh"
source "$HOME/lab_utils/bash/functions/bmc_utils.sh"
source "$HOME/lab_utils/bash/functions/fastq_utils.sh"

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

    # Ensure log file is writable
    if ! touch "$log_file" 2>/dev/null; then
        echo "ERROR: Unable to create or write to log file: $log_file. Check permissions or quotas."
        return 1
    fi

    echo "$log_file"
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
    if ! move_fastq_files_to_current_directory "${paths#*:}" "$log_file"; then
        return 1
    fi

    # Cleanup
    if ! clean_experiment_directory "${paths#*:}" "$log_file"; then
        return 1
    fi

    log_info "Download process completed successfully" "$log_file"
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
