#!/bin/bash
# bash/scripts/transfer_bmc_data.sh

# Configuration
REMOTE_HOST="luria.mit.edu"
REMOTE_USER="luised94"
REMOTE_PATH="~/data"
LOG_DIR="$HOME/logs"

# Ensure log directory exists
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/transfer_$(date +%Y%m%d_%H%M%S).log"

# Logging function
log_message() {
    local level="$1"
    local message="$2"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$level] $message" | tee -a "$LOG_FILE"
}

# Check network connectivity
check_connection() {
    if ! ping -c 1 "$REMOTE_HOST" &> /dev/null; then
        log_message "ERROR" "Cannot connect to $REMOTE_HOST"
        log_message "INFO" "Please ensure:"
        log_message "INFO" "1. VPN is connected"
        log_message "INFO" "2. You can connect via: ssh -A -Y $REMOTE_USER@$REMOTE_HOST"
        return 1
    fi
    return 0
}

# Validate local directory
validate_directory() {
    local dir="$1"
    if [ ! -d "$dir" ]; then
        log_message "ERROR" "Directory not found: $dir"
        return 1
    }
    
    # Check for required subdirectories
    for subdir in "documentation" "fastq" "logs"; do
        if [ ! -d "$dir/$subdir" ]; then
            log_message "ERROR" "Required subdirectory missing: $subdir"
            return 1
        fi
    done
    return 0
}

# Transfer data
transfer_data() {
    local source_dir="$1"
    local experiment_id=$(basename "$source_dir")
    
    log_message "INFO" "Starting transfer of $experiment_id"
    
    rsync -avzP --stats \
        "$source_dir/" \
        "$REMOTE_USER@$REMOTE_HOST:$REMOTE_PATH/$experiment_id/" \
        2>&1 | tee -a "$LOG_FILE"
    
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        log_message "INFO" "Transfer completed successfully"
        return 0
    else
        log_message "ERROR" "Transfer failed"
        return 1
    fi
}

# Verify transfer
verify_transfer() {
    local source_dir="$1"
    local experiment_id=$(basename "$source_dir")
    
    log_message "INFO" "Verifying transfer"
    
    # Compare file counts
    local local_count=$(find "$source_dir" -type f | wc -l)
    local remote_count=$(ssh "$REMOTE_USER@$REMOTE_HOST" "find $REMOTE_PATH/$experiment_id -type f | wc -l")
    
    if [ "$local_count" -eq "$remote_count" ]; then
        log_message "INFO" "Verification successful: $local_count files transferred"
        return 0
    else
        log_message "ERROR" "Verification failed: Local=$local_count Remote=$remote_count"
        return 1
    fi
}

# Main function
main() {
    if [ $# -ne 1 ]; then
        log_message "ERROR" "Usage: $0 <experiment_directory>"
        exit 1
    fi
    
    local source_dir="$1"
    
    # Run checks
    check_connection || exit 1
    validate_directory "$source_dir" || exit 1
    
    # Transfer and verify
    transfer_data "$source_dir" || exit 1
    verify_transfer "$source_dir" || exit 1
    
    log_message "INFO" "Next steps:"
    log_message "INFO" "1. Login to cluster: ssh -A -Y $REMOTE_USER@$REMOTE_HOST"
    log_message "INFO" "2. Run: bash ~/lab_utils/bash/scripts/download_bmc_data.sh"
}

# Execute main function
main "$@"
