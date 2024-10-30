#!/bin/bash
# bash/scripts/transfer_bmc_experiment_to_luria.sh
# Requires password twice, for the transfer and verification.
# Run after setup_bmc_experiment.R

# Source dependencies
source "$HOME/lab_utils/bash/functions/logging_utils.sh"

#' Check Network Connectivity
#' @param remote_host Character Remote host address
#' @return Integer 0 if successful, 1 otherwise
check_connection() {
    local remote_host="$1"
    
    if ! ping -c 1 "$remote_host" &> /dev/null; then
        log_error "Cannot connect to $remote_host"
        log_info "Please ensure:"
        log_info "1. VPN is connected"
        log_info "2. You can connect via: ssh -A -Y ${PROJECT_CONFIG[REMOTE_USER]}@$remote_host"
        return 1
    fi
    return 0
}

#' Validate Local Directory Structure
#' @param dir Character Directory to validate
#' @return Integer 0 if valid, 1 otherwise
validate_directory() {
    local dir="$1"
    local log_file="$2"
    
    if [ ! -d "$dir" ]; then
        log_error "Directory not found: $dir" "${log_file}"
        log_info "Provide full path." "${log_file}"
        return 1
    fi
    
    # Check required subdirectories
    for subdir in ${PROJECT_CONFIG[REQUIRED_DIRS]}; do
        if [ ! -d "$dir/$subdir" ]; then
            log_error "Required subdirectory missing: $subdir"
            return 1
        fi
    done
    return 0
}

#' Transfer Data to Remote Host
#' @param source_dir Character Source directory
#' @param log_file Character Log file path
#' @return Integer 0 if successful, 1 otherwise
transfer_data() {
    local source_dir="$1"
    local log_file="$2"
    local experiment_id=$(basename "$source_dir")
    
    log_info "Starting transfer of $experiment_id" "$log_file"
    
    rsync -avzP --stats \
        "$source_dir/" \
        "${PROJECT_CONFIG[REMOTE_USER]}@${PROJECT_CONFIG[REMOTE_HOST]}:${PROJECT_CONFIG[REMOTE_PATH]}/$experiment_id/" \
        2>&1 | tee -a "$log_file"
    
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        log_info "Transfer completed successfully" "$log_file"
        return 0
    else
        log_error "Transfer failed" "$log_file"
        return 1
    fi
}

#' Verify Transfer Completion
#' @param source_dir Character Source directory
#' @param log_file Character Log file path
#' @return Integer 0 if verified, 1 otherwise
verify_transfer() {
    local source_dir="$1"
    local log_file="$2"
    local experiment_id=$(basename "$source_dir")
    
    log_info "Verifying transfer" "$log_file"
    
    local local_count=$(find "$source_dir" -type f | wc -l)
    local remote_count=$(ssh "${PROJECT_CONFIG[REMOTE_USER]}@${PROJECT_CONFIG[REMOTE_HOST]}" \
        "find ${PROJECT_CONFIG[REMOTE_PATH]}/$experiment_id -type f | wc -l")
    
    if [ "$local_count" -eq "$remote_count" ]; then
        log_info "Verification successful: $local_count files transferred" "$log_file"
        return 0
    else
        log_error "Verification failed: Local=$local_count Remote=$remote_count" "$log_file"
        return 1
    fi
}

#' Main Function
#' @param args Array Script arguments
#' @return None
transfer_bmc_experiment_to_luria_main() {
    # Initialize logging
    local log_file
    log_file="$(initialize_logging "transfer_bmc_experiment")"
    
    if [ $# -ne 1 ]; then
        log_error "Usage: $0 <experiment_directory>" "$log_file"
        exit 1
    fi
    
    local source_dir="$1"
    
    # Run checks
    check_connection "${PROJECT_CONFIG[REMOTE_HOST]}" || exit 1
    validate_directory "$source_dir" "$log_file" || exit 1
    
    # Transfer and verify
    transfer_data "$source_dir" "$log_file" || exit 1
    verify_transfer "$source_dir" "$log_file" || exit 1
    
    log_info "Next steps:" "$log_file"
    log_info "1. Login to cluster: ssh -A -Y ${PROJECT_CONFIG[REMOTE_USER]}@${PROJECT_CONFIG[REMOTE_HOST]}" "$log_file"
    log_info "2. Run: bash ~/lab_utils/bash/scripts/download_bmc_data.sh" "$log_file"
}

# Execute if run as script
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    transfer_bmc_experiment_to_luria_main "$@"
fi
