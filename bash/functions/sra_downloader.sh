#!/bin/bash

source "../config/sra_config.sh"

function validate_input() {
    local download_dir="$1"
    
    log_info "Validating input parameters"
    
    if [ -z "$download_dir" ]; then
        log_error "Download directory not specified"
        return 1
    }
    
    local full_path="${SRA_CONFIG[DATA_DIR]}/$download_dir"
    
    if [ ! -d "$full_path" ]; then
        log_info "Creating directory: $full_path"
        mkdir -p "$full_path" || {
            log_error "Failed to create directory: $full_path"
            return 1
        }
    }
    
    echo "$full_path"
}

function construct_download_url() {
    local accession="$1"
    local base_url="${SRA_CONFIG[BASE_URL]}"
    
    echo "${base_url}${accession:0:6}/${accession}/${accession}.fastq.gz"
}

function verify_url() {
    local url="$1"
    
    log_info "Verifying URL: $url"
    
    if ! curl --head --silent --fail "$url" >/dev/null; then
        log_error "URL not accessible: $url"
        return 1
    }
    
    return 0
}

function download_file() {
    local url="$1"
    local output_file="$2"
    
    log_info "Downloading: $url -> $output_file"
    
    if ! wget --quiet --show-progress --output-document="$output_file" "$url"; then
        log_error "Download failed: $url"
        return 1
    }
    
    log_info "Download complete: $output_file"
    return 0
}

function concatenate_files() {
    local output_dir="$1"
    local output_file="$2"
    local files=("${@:3}")
    
    log_info "Concatenating files to: $output_file"
    
    for file in "${files[@]}"; do
        if [ ! -f "$output_dir/$file" ]; then
            log_error "File not found: $file"
            return 1
        }
        cat "$output_dir/$file" >> "$output_file" || {
            log_error "Failed to concatenate: $file"
            return 1
        }
    done
    
    log_info "Concatenation complete"
    return 0
}

function decompress_file() {
    local file="$1"
    
    log_info "Decompressing: $file"
    
    if ! gunzip "$file"; then
        log_error "Decompression failed: $file"
        return 1
    }
    
    log_info "Decompression complete"
    return 0
}
