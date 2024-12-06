#!/bin/bash

source "../functions/ngs_file_manager.sh"

function setup_qc_directory() {
    local base_dir="$1"
    local qc_dir="${base_dir}/${QC_CONFIG[OUTPUT_DIR]}"
    
    log_info "Setting up QC directory: $qc_dir"
    mkdir -p "$qc_dir" || {
        log_error "Failed to create QC directory: $qc_dir"
        return 1
    }
    
    echo "$qc_dir"
}

function run_fastqc() {
    local file_path="$1"
    local output_dir="$2"
    
    if [ ! -f "$file_path" ]; then
        log_error "File not found: $file_path"
        return 1
    fi
    
    log_info "Running FASTQC on: $file_path"
    if ! fastqc --outdir="$output_dir" "$file_path"; then
        log_error "FASTQC failed for: $file_path"
        return 1
    fi
    
    log_info "FASTQC completed for: $file_path"
    return 0
}
