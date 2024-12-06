#!/bin/bash
# functions moved

set -euo pipefail

source "../functions/archive_handler.sh"

function show_usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS] <directory>
Unzips FASTQC result files in specified directory

Options:
    -f, --force       Skip confirmation
    -k, --keep        Keep ZIP files after extraction
    -b, --batch SIZE  Process in batches of SIZE files
    -h, --help        Show this help message

Example:
    $(basename "$0") experiment_20240101
EOF
}

function confirm_operation() {
    local files="$1"
    local timeout="${OPERATION_DEFAULTS[CONFIRM_TIMEOUT]}"
    
    echo "Files to process:"
    echo "$files"
    echo
    
    read -t "$timeout" -p "Proceed with unzipping these files? (y/n): " -r || {
        echo
        log_error "Confirmation timed out"
        return 1
    }
    
    [[ $REPLY =~ ^[Yy]$ ]]
}

function validate_directory() {
    local dir_name="$1"
    local base_dir="${QC_PATHS[BASE_DIR]}"
    
    log_info "Validating directory: $dir_name"
    
    local full_path=$(find -H "$base_dir" \
                      -maxdepth "${QC_PATHS[MAX_DEPTH]}" \
                      -type d \
                      -name "$dir_name")
    
    if [ -z "$full_path" ]; then
        log_error "Directory not found: $dir_name"
        return 1
    }
    
    echo "$full_path"
}

function main() {
    local force=false
    local preserve=${OPERATION_DEFAULTS[PRESERVE_ZIP]}
    local batch_size=${OPERATION_DEFAULTS[UNZIP_BATCH_SIZE]}
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -f|--force) force=true; shift ;;
            -k|--keep) preserve=true; shift ;;
            -b|--batch) batch_size="$2"; shift 2 ;;
            -h|--help) show_usage; exit 0 ;;
            *) break ;;
        esac
    done
    
    if [ $# -ne 1 ]; then
        log_error "Directory name required"
        show_usage
        exit 1
    fi
    
    local exp_dir=$(validate_directory "$1") || exit 1
    local qc_dir="$exp_dir/${QC_PATHS[QC_SUBDIR]}"
    
    mapfile -t zip_files < <(find_zip_files "$qc_dir") || {
        log_error "No ZIP files found in: $qc_dir"
        exit 1
    }
    
    if [ "$force" = false ]; then
        confirm_operation "${zip_files[*]}" || exit 1
    fi
    
    process_zip_files "${zip_files[@]}"
}

main "$@"
