#!/bin/bash
# functions moved

set -euo pipefail

source "../functions/sra_downloader.sh"

function show_usage() {
    cat << EOF
Usage: $(basename "$0") <directory>
Downloads Eaton 2010 paper data for control and comparison

Arguments:
    directory    Target directory name (relative to ${SRA_CONFIG[DATA_DIR]})

Study Information:
    BioProject: ${EATON_2010[BIOPROJECT]}
    Description: ${EATON_2010[DESCRIPTION]}
EOF
}

function main() {
    if [ $# -ne 1 ]; then
        show_usage
        exit 1
    }
    
    local output_dir=$(validate_input "$1") || exit 1
    local combined_output="${output_dir}/${EATON_2010[OUTPUT_FILE]}.gz"
    local downloaded_files=()
    
    for file_name in "${!SAMPLES[@]}"; do
        local accession="${SAMPLES[$file_name]}"
        local url=$(construct_download_url "$accession")
        local output_file="${output_dir}/$file_name"
        
        verify_url "$url" || continue
        download_file "$url" "$output_file" || continue
        
        downloaded_files+=("$file_name")
    done
    
    if [ ${#downloaded_files[@]} -gt 0 ]; then
        concatenate_files "$output_dir" "$combined_output" "${downloaded_files[@]}" || exit 1
        decompress_file "$combined_output" || exit 1
    else
        log_error "No files were downloaded successfully"
        exit 1
    fi
    
    log_info "Processing complete"
}

main "$@"
