#!/bin/bash

source "../config/file_management_config.sh"

function build_file_patterns() {
    local patterns=()
    for type in "${!FILE_TYPES[@]}"; do
        for ext in ${FILE_TYPES[$type]}; do
            patterns+=("-o" "-name" "*.$ext")
        done
    done
    echo "${patterns[@]:1}" # Remove first -o
}

function check_disk_space() {
    local target_dir="$1"
    local min_space="${2:-${DEFAULTS[MIN_FREE_SPACE_GB]}}"
    
    local available_space=$(df -BG "$target_dir" | awk 'NR==2 {print $4}' | sed 's/G//')
    if (( available_space < min_space )); then
        log_error "Insufficient disk space. Available: ${available_space}GB, Required: ${min_space}GB"
        return 1
    }
    return 0
}

function validate_directories() {
    local target_dir="$1"
    local search_dir="$2"
    
    if [[ ! -d "$search_dir" ]]; then
        log_error "Search directory does not exist: $search_dir"
        return 1
    }
    
    if [[ ! -d "$target_dir" ]]; then
        log_info "Creating target directory: $target_dir"
        mkdir -p "$target_dir"
    }
    
    if [[ ! -w "$target_dir" ]]; then
        log_error "Target directory not writable: $target_dir"
        return 1
    }
}

function move_ngs_files() {
    local target_dir="$1"
    local search_dir="${2:-.}"
    local max_depth="${3:-${DEFAULTS[MAX_DEPTH]}}"
    local batch_size="${4:-${DEFAULTS[BATCH_SIZE]}}"
    
    log_info "Starting NGS file movement operation"
    
    validate_directories "$target_dir" "$search_dir" || return 1
    check_disk_space "$target_dir" || return 1
    
    local file_patterns=($(build_file_patterns))
    
    find "$search_dir" \
        -maxdepth "$max_depth" \
        -type f \
        \( "${file_patterns[@]}" \) \
        \( -path "*/code/*" -o -path "*/script*/*" \) \
        -print0 | 
    while IFS= read -r -d '' file; do
        local file_type=$(determine_file_type "$file")
        local target_subdir="$target_dir/$file_type"
        
        mkdir -p "$target_subdir"
        log_info "Moving $file to $target_subdir"
        mv "$file" "$target_subdir/" || log_error "Failed to move: $file"
    done
}

function determine_file_type() {
    local file="$1"
    local extension="${file##*.}"
    
    for type in "${!FILE_TYPES[@]}"; do
        if [[ " ${FILE_TYPES[$type]} " =~ " $extension " ]]; then
            echo "$type"
            return
        }
    done
    echo "OTHER"
}

function analyze_file_distribution() {
    local search_dir="${1:-.}"
    
    log_info "Analyzing file distribution in $search_dir"
    
    find "$search_dir" \
        -type d \( -path '*/lib/*' -o -path '*/renv/*' -o -path '*/git/*' \) -prune \
        -o -type f -exec du -a {} + | 
        sort -nr | 
        head -n 1000 | 
        awk -F/ '{print $NF}' | 
        rev | cut -d. -f1 | rev | 
        sort | uniq -c | 
        sort -nr
}

# Add to existing file management functions
function find_fastq_files() {
    local base_dir="$1"
    local exclude_pattern="${QC_CONFIG[EXCLUDE_PATTERNS]}"
    
    log_info "Finding FASTQ files in: $base_dir"
    find "$base_dir" -type f -name "*.fastq" ! \( -name "*unmapped*" -o -name "processed_*" \)
}

function find_processed_fastq() {
    local base_dir="$1"
    local pattern="${QC_CONFIG[PROCESSED_PATTERN]}"
    
    log_info "Finding processed FASTQ files in: $base_dir"
    find "$base_dir" -type f -name "$pattern"
}
