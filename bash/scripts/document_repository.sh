#!/bin/bash

set -euo pipefail

source "../functions/file_processor.sh"

function show_usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]
Generates XML documentation of repository contents

Options:
    -q, --quiet       Suppress verbose output
    -o, --output FILE Output file name (default: ${DOC_CONFIG[OUTPUT_FILE]})
    -h, --help        Show this help message
EOF
}

function main() {
    local verbose=${DOC_CONFIG[DEFAULT_VERBOSE]}
    local output_file=${DOC_CONFIG[OUTPUT_FILE]}
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -q|--quiet) verbose=0; shift ;;
            -o|--output) output_file="$2"; shift 2 ;;
            -h|--help) show_usage; exit 0 ;;
            *) log_error "Unknown option: $1"; show_usage; exit 1 ;;
        esac
    done
    
    local root_dir=$(get_repository_root)
    
    create_xml_header "$output_file"
    
    for category in "${!FILE_TYPES[@]}"; do
        for type in ${FILE_TYPES[$category]}; do
            process_file_type "$type" "$output_file" "$root_dir"
        done
    done
    
    close_xml_document "$output_file"
    
    log_info "Documentation generated: $output_file"
}

main "$@"
