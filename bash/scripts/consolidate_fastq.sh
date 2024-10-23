#!/bin/bash
# functions moved

set -eo pipefail

source "../functions/fastq_processor.sh"

function show_usage() {
    cat << EOF
Usage: $(basename "$0") <experiment_name>
Consolidates multiple FASTQ files per sample into single files.

Arguments:
    experiment_name    Name of the experiment (without trailing slash)

Example: 
    $(basename "$0") "240808Bel"
EOF
}

function main() {
    if [ $# -ne 1 ]; then
        show_usage
        exit 1
    }
    
    local experiment_dir=$(validate_input "$1") || exit 1
    local output_dir=$(setup_output_directory "$experiment_dir") || exit 1
    
    mapfile -t fastq_files < <(find_fastq_files "$experiment_dir")
    local initial_count=${#fastq_files[@]}
    
    mapfile -t unique_ids < <(extract_unique_ids fastq_files)
    
    for id in "${unique_ids[@]}"; do
        process_fastq_file "$id" "$output_dir" fastq_files || exit 1
    done
    
    verify_processing "$experiment_dir" "$1" "$initial_count" "${#unique_ids[@]}"
}

main "$@"
