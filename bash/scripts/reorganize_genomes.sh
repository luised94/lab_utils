#!/bin/bash
# functions moved

set -euo pipefail

source "../functions/genome_organizer.sh"

function show_usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]
Reorganizes reference genome directories

Options:
    -d, --dir DIR     Process specific directory
    -a, --all         Process all genome directories
    -s, --standardize Standardize chromosome names
    -h, --help        Show this help message
EOF
}

function process_genome_directory() {
    local dir="$1"
    
    log_info "Processing directory: $dir"
    
    local assembly_report="${dir}/${GENOME_PATHS[NCBI_DATA]}/${GENOME_PATHS[ASSEMBLY_REPORT]}"
    
    local organism_name=$(extract_organism_name "$assembly_report") || return 1
    
    if ! reorganize_genome_files "$dir" "$organism_name"; then
        log_error "Failed to reorganize: $dir"
        return 1
    }
    
    if [ -d "$dir" ] && [ "$dir" != "$organism_name" ]; then
        mv "$dir" "$organism_name" || {
            log_error "Failed to rename directory to: $organism_name"
            return 1
        }
    }
    
    log_info "Successfully processed: $organism_name"
    return 0
}

function main() {
    local process_all=false
    local standardize=false
    local target_dir=""
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -d|--dir) target_dir="$2"; shift 2 ;;
            -a|--all) process_all=true; shift ;;
            -s|--standardize) standardize=true; shift ;;
            -h|--help) show_usage; exit 0 ;;
            *) log_error "Unknown option: $1"; show_usage; exit 1 ;;
        esac
    done
    
    if $process_all; then
        mapfile -t dirs < <(find . -maxdepth 1 -type d -name "GC[AF]_*")
    elif [ -n "$target_dir" ]; then
        dirs=("$target_dir")
    else
        log_error "No directory specified"
        show_usage
        exit 1
    fi
    
    for dir in "${dirs[@]}"; do
        process_genome_directory "$dir" || continue
    done
    
    log_info "Genome reorganization completed"
}

main "$@"
