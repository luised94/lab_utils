#!/bin/bash
# functions moved

set -euo pipefail

source "../functions/bmc_transfer.sh"
source "../functions/file_organizer.sh"

function show_usage() {
    cat << EOF
Usage: $(basename "$0") <bmc_server> <experiment_directory>
Downloads FASTQ data from BMC servers

Arguments:
    bmc_server           BMC server name (e.g., bmc-pub17)
    experiment_directory Experiment directory name

Example:
    $(basename "$0") bmc-pub17 240808Bel
EOF
}

function cleanup_data() {
    log_info "Starting cleanup process"
    
    for pattern in "${CLEANUP_PATTERNS[DIRS][@]}"; do
        safe_remove "dir" "$pattern"
    done
    
    for pattern in "${CLEANUP_PATTERNS[FILES][@]}"; do
        safe_remove "file" "$pattern"
    done
}

function main() {
    if [ $# -ne 2 ]; then
        show_usage
        exit 1
    }
    
    local server="$1"
    local experiment="$2"
    
    if ! download_bmc_data "$server" "$experiment"; then
        log_error "Download failed"
        exit 1
    }
    
    local target_dir="${BMC_CONFIG[LOCAL_BASE]}/$experiment/${BMC_CONFIG[FASTQ_DIR]}"
    
    if ! organize_fastq_files "$target_dir"; then
        log_error "File organization failed"
        exit 1
    }
    
    cleanup_data
    
    log_info "Process completed successfully"
    log_info "Verify files with: find $target_dir -type f -name '*.fastq' | wc -l"
}

main "$@"
