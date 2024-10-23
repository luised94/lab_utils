#!/bin/bash

set -euo pipefail

source "../functions/feature_data_handler.sh"

function show_usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]
Downloads feature data from Rossi 2021 paper

Options:
    -d, --dir DIR     Specify download directory
                      (default: ${FEATURE_DATA[BASE_DIR]})
    -h, --help        Show this help message

Description:
    Downloads and sets up feature data from the Rossi 2021 paper
    for categorical analysis and plot tracking.
EOF
}

function main() {
    local download_dir="${FEATURE_DATA[BASE_DIR]}"
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -d|--dir) download_dir="$2"; shift 2 ;;
            -h|--help) show_usage; exit 0 ;;
            *) log_error "Unknown option: $1"; show_usage; exit 1 ;;
        esac
    done
    
    validate_git || exit 1
    
    local base_dir=$(setup_feature_directory) || exit 1
    
    if ! download_rossi_data "$base_dir"; then
        log_error "Failed to download Rossi data"
        exit 1
    fi
    
    log_info "Feature data download completed successfully"
}

main "$@"
