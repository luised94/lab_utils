#!/bin/bash
# bash/scripts/fastq/download_bmc_data.sh

#' Download BMC Data Script
#' @description Download and process BMC FASTQ data

echo "Initializing BMC data download"
echo "Script location: ${BASH_SOURCE[0]}"
echo "Working directory: $(pwd)"

# Find repository root using git
repo_root=$(git rev-parse --show-toplevel 2>/dev/null) || {
    echo "? Not in a git repository"
    exit 1
}

# Initialize environment
source "$repo_root/bash/core/initialize_lab_environment.sh" || {
    echo "? Failed to initialize environment"
    exit 1
}


# Load required modules
echo "Loading required modules"
for module in "fastq/bmc_handler" "fastq/fastq_processor"; do
    echo "Loading: $module"
    if ! load_lab_module "$module"; then
        log_error "Failed to load module: $module"
        exit 1
    fi
    echo "Loaded successfully"
done
#' Download BMC Data Main Function
#' @param experiment_id Character Experiment identifier
#' @return Integer 0 if successful
download_bmc_data_main() {
    local experiment_id="$1"
    local log_file
    
    # Initialize logging
    log_file=$(initialize_logging "download_bmc_data") || {
        echo "? Failed to initialize logging"
        return 1
    }
    
    log_info "Starting BMC data download" "$log_file"
    log_debug "Environment verification" "$log_file"
    log_debug "LAB_UTILS_ROOT: $LAB_UTILS_ROOT" "$log_file"
    log_debug "Experiment ID: $experiment_id" "$log_file"
    log_debug "Log file: $log_file" "$log_file"
    
    # Verify host
    if ! verify_host; then
        log_error "Host verification failed" "$log_file"
        return 1
    fi
    
    # Validate paths
    local paths
    log_debug "Validating BMC paths" "$log_file"
    if ! paths=$(validate_bmc_paths "$experiment_id" "$log_file"); then
        return 1
    fi
    
    local bmc_path=${paths%:*}
    local local_path=${paths#*:}

    log_debug " BMC path: $bmc_path" "$log_file"
    log_debug " Local path: $local_path" "$log_file"
    
    # Download data
    if ! download_from_bmc "$paths" "$log_file"; then
        return 1
    fi

    # Organize files
    if ! move_fastq_files_to_current_directory "${local_path}" "$log_file"; then
        return 1
    fi

    # Cleanup
    if ! clean_experiment_directory "${local_path}" "$log_file"; then
        return 1
    fi

    log_info "Download process completed successfully" "$log_file"
    return 0
}

# Show usage information
show_usage() {
    cat << EOF
Usage: $(basename "$0") <experiment_id>
Arguments:
    experiment_id    Experiment identifier (format: YYMMDD'Bel')
Example:
    $(basename "$0") 241010Bel
EOF
}

# Execute if run as script
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    download_bmc_data_main "$@"
fi
