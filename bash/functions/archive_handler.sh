#!/bin/bash

source "../config/quality_control_config.sh"

function find_zip_files() {
    local base_dir="$1"
    local pattern="${FILE_PATTERNS[FASTQC_ZIP]}"
    
    log_info "Searching for ZIP files in: $base_dir"
    
    local files=$(find "$base_dir" -type f -name "$pattern")
    
    if [ -z "$files" ]; then
        log_warning "No ZIP files found"
        return 1
    }
    
    echo "$files"
}

function verify_zip_file() {
    local zip_file="$1"
    
    if ! unzip -t "$zip_file" >/dev/null 2>&1; then
        log_error "Invalid or corrupted ZIP file: $zip_file"
        return 1
    }
    
    return 0
}

function unzip_in_place() {
    local zip_file="$1"
    local preserve="${2:-${OPERATION_DEFAULTS[PRESERVE_ZIP]}}"
    
    local dir=$(dirname "$zip_file")
    local filename=$(basename "$zip_file")
    
    log_info "Unzipping: $filename"
    
    (
        cd "$dir" || {
            log_error "Failed to change to directory: $dir"
            return 1
        }
        
        if ! unzip -o "$filename"; then
            log_error "Failed to unzip: $filename"
            return 1
        }
        
        if [ "$preserve" = false ]; then
            rm "$filename" || log_warning "Failed to remove ZIP file: $filename"
        fi
    )
}

function process_zip_files() {
    local files=("$@")
    local total=${#files[@]}
    local success=0
    local failed=0
    
    log_info "Processing $total ZIP files"
    
    for file in "${files[@]}"; do
        if verify_zip_file "$file"; then
            if unzip_in_place "$file"; then
                ((success++))
            else
                ((failed++))
            fi
        else
            ((failed++))
        fi
    done
    
    log_info "Processed files - Success: $success, Failed: $failed"
    return $((failed > 0))
}
