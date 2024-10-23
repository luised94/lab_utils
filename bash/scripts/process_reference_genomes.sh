#!/bin/bash

set -euo pipefail

source "../functions/ncbi_handler.sh"
source "../functions/genome_indexer.sh"  # From previous implementation

function show_usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]
Downloads and processes reference genomes from NCBI

Options:
    -d, --download    Download genomes only
    -i, --index       Build indices only
    -a, --all         Download and index genomes
    -h, --help        Show this help message
EOF
}

function main() {
    local do_download=false
    local do_index=false
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -d|--download) do_download=true; shift ;;
            -i|--index) do_index=true; shift ;;
            -a|--all) do_download=true; do_index=true; shift ;;
            -h|--help) show_usage; exit 0 ;;
            *) log_error "Unknown option: $1"; show_usage; exit 1 ;;
        esac
    done
    
    if ! $do_download && ! $do_index; then
        log_error "No operation specified"
        show_usage
        exit 1
    }
    
    validate_ncbi_tools || exit 1
    
    if $do_download; then
        log_info "Starting genome downloads"
        process_genome_batch "${NCBI_CONFIG[ACCESSIONS][@]}" || exit 1
    fi
    
    if $do_index; then
        log_info "Starting genome indexing"
        build_genome_indices || exit 1  # From previous implementation
    fi
    
    log_info "All operations completed successfully"
}

main "$@"
