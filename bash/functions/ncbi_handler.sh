#!/bin/bash

source "../config/genome_config.sh"

function validate_ncbi_tools() {
    if ! command -v datasets &>/dev/null; then
        log_error "NCBI datasets tool not found"
        return 1
    }
    return 0
}

function download_genome() {
    local accession="$1"
    local base_dir="${GENOME_CONFIG[BASE_DIR]}"
    local timestamp=$(date +"%Y%m%d_%H%M%S")
    local target_dir="${base_dir}/${accession}_${timestamp}"
    
    log_info "Downloading genome: $accession"
    
    mkdir -p "$target_dir" || {
        log_error "Failed to create directory: $target_dir"
        return 1
    }
    
    if ! datasets download genome accession "$accession" \
        --include "${NCBI_CONFIG[DOWNLOAD_INCLUDES]}" \
        --filename "${target_dir}/${accession}.zip"; then
        log_error "Download failed for accession: $accession"
        return 1
    }
    
    log_info "Extracting files for: $accession"
    if ! unzip -q "${target_dir}/${accession}.zip" -d "$target_dir"; then
        log_error "Extraction failed for: $accession"
        return 1
    }
    
    rm "${target_dir}/${accession}.zip"
    log_info "Successfully processed: $accession"
    return 0
}

function process_genome_batch() {
    local -a accessions=("$@")
    local success_count=0
    local fail_count=0
    
    for accession in "${accessions[@]}"; do
        if download_genome "$accession"; then
            ((success_count++))
        else
            ((fail_count++))
            log_warning "Failed to process: $accession"
        fi
    done
    
    log_info "Processing complete. Success: $success_count, Failed: $fail_count"
    return $((fail_count > 0))
}
