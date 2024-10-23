#!/bin/bash

source "../config/slurm_config.sh"
# Assuming logging utilities are already sourced

function find_slurm_files() {
    local search_dir="${1:-.}"
    local pattern="${2:-$SLURM_OUTPUT_PATTERN}"
    local max_depth="${3:-$MAX_DEPTH_SEARCH}"

    log_info "Searching for SLURM files in: $search_dir"
    
    if [[ ! -d "$search_dir" ]]; then
        log_error "Directory does not exist: $search_dir"
        return 1
    }

    find "$search_dir" \
        -maxdepth "$max_depth" \
        -type f \
        -name "$pattern" \
        -printf "%T@ %p\n" | \
        sort -nr | \
        cut -d' ' -f2-
}

function organize_slurm_files() {
    local source_dir="${1:-.}"
    local target_dir="${2:-$DEFAULT_SLURM_LOG_DIR}"
    
    log_info "Organizing SLURM files from $source_dir to $target_dir"
    
    # Create target directory if it doesn't exist
    mkdir -p "$target_dir"
    
    # Move files to organized structure
    while IFS= read -r file; do
        local job_id=$(basename "$file" | sed 's/slurm-\([0-9]*\).out/\1/')
        local date_str=$(stat -c %y "$file" | cut -d' ' -f1)
        local year_month=$(date -d "$date_str" +%Y/%m)
        local target_subdir="$target_dir/$year_month"
        
        mkdir -p "$target_subdir"
        mv "$file" "$target_subdir/"
        log_info "Moved $file to $target_subdir/"
    done < <(find_slurm_files "$source_dir")
}

function cleanup_old_logs() {
    local log_dir="${1:-$DEFAULT_SLURM_LOG_DIR}"
    local max_age="${2:-$MAX_LOG_AGE_DAYS}"
    
    log_info "Cleaning up SLURM logs older than $max_age days"
    
    find "$log_dir" \
        -type f \
        -name "$SLURM_OUTPUT_PATTERN" \
        -mtime +"$max_age" \
        -exec rm {} \;
}

function check_large_files() {
    local search_dir="${1:-.}"
    local max_size="${2:-$MAX_FILE_SIZE_MB}"
    
    log_info "Checking for large SLURM output files"
    
    find "$search_dir" \
        -type f \
        -name "$SLURM_OUTPUT_PATTERN" \
        -size +"${max_size}M" \
        -exec ls -lh {} \; | \
    while read -r line; do
        log_warning "Large SLURM output file: $line"
    done
}
