#!/bin/bash
# bash/functions/fastq_processor.sh

source "$HOME/lab_utils/bash/functions/logging_utils.sh"

#' Validate FASTQ Processing Input
#' @param experiment_id Character Experiment identifier
#' @param log_file Character Log file path
#' @return String Experiment directory path
validate_fastq_input() {
    local experiment_id="$1"
    local log_file="$2"
    {
        echo "DEBUG:"
        echo "$experiment_id"
        echo "$log_file"
    } >&2
    log_info "Validating input for experiment: $experiment_id" "$log_file" >&2
    
    # Validate experiment ID format
    if [[ ! "$experiment_id" =~ ^[0-9]{6}Bel$ ]]; then
        log_error "Invalid experiment ID format: $experiment_id" "$log_file"
        return 1
    fi
    
    local experiment_dir="$HOME/data/$experiment_id"
    if [[ ! -d "$experiment_dir" ]]; then
        log_error "Experiment directory not found: $experiment_dir" "$log_file"
        return 1
    fi
    
    echo -n "$experiment_dir"
}

#' Setup FASTQ Processing Directories
#' @param experiment_dir Character Experiment directory
#' @param log_file Character Log file path
#' @return String Output directory path
setup_fastq_directories() {
    local experiment_dir="$1"
    local log_file="$2"
    
{
    echo "DEBUG:"
    echo "$experiment_dir"
    echo "$log_file"
} >&2
    local output_dir="${experiment_dir}/fastq"
    log_trace "$output_dir" >&2
    log_info "Setting up output directory: $output_dir" "$log_file" >&2
    
    mkdir -p "$output_dir" || {
        log_error "Failed to create output directory: $output_dir" "$log_file"
        return 1
    }
    
    echo -n "$output_dir"
}

#' Find FASTQ Files
#' @param experiment_dir Character Experiment directory
#' @param log_file Character Log file path
#' @return Array FASTQ file paths

find_fastq_files() {
    local experiment_dir="$1"
    local log_file="$2"
    
    log_info "Searching for FASTQ files in: $experiment_dir/${PROJECT_CONFIG[FASTQ_DIR]}" "$log_file" >&2
    
    # Validate directory exists
    if [[ ! -d "$experiment_dir/${PROJECT_CONFIG[FASTQ_DIR]}" ]]; then
        log_error "FASTQ directory not found: $experiment_dir/${PROJECT_CONFIG[FASTQ_DIR]}" "$log_file"
        return 1
    fi
    
    # Build exclude patterns array
    local -a exclude_patterns=()
    for pattern in ${PROJECT_CONFIG[FASTQ_EXCLUDE]}; do
        exclude_patterns+=( "!" "-name" "$pattern" )
    done
    
    # Execute find with error handling
    local find_output
    if ! find_output=$(find "$experiment_dir/${PROJECT_CONFIG[FASTQ_DIR]}" -type f \
        -name "${PROJECT_CONFIG[FASTQ_PATTERN]}" \
        "${exclude_patterns[@]}" 2>&1 | sort); then
        log_error "Find command failed: $find_output" "$log_file"
        return 1
    fi
    
    # Check if any files were found
    if [[ -z "$find_output" ]]; then
        log_warning "No FASTQ files found matching pattern: ${PROJECT_CONFIG[FASTQ_PATTERN]}" "$log_file"
        return 0
    fi
    
    echo -n "$find_output"
}

#' Process FASTQ Files
#' @param experiment_dir Character Experiment directory
#' @param output_dir Character Output directory
#' @param log_file Character Log file path
#' @return Integer 0 if successful
process_fastq_files() {
    local experiment_dir="$1"
    local output_dir="$2"
    local log_file="$3"
{
    echo "DEBUG:"
    echo "$experiment_dir"
    echo "$log_file"
    echo "$output_dir"
} >&2
    
    local -a fastq_files
    mapfile -t fastq_files < <(find_fastq_files "$experiment_dir" "$log_file")
    
    local initial_count=${#fastq_files[@]}
    log_info "Found $initial_count FASTQ files to process" "$log_file" >&2

    local debug_file=${fastq_files[1]}
    local debug_basename=$(basename debug_file)
    local debug_id=$(echo "$debug_basename" | grep -oP "${PROJECT_CONFIG[BMC_FASTQ_ID_PATTERN]}\K[^_]*")
    local debug_output_file="$output_dir/${id}${PROJECT_CONFIG[FASTQ_SUFFIX]}"
{
    echo "DEBUG:"
    echo "$debug_file"
    echo "$debug_basename"
    echo "$debug_id"
    echo "$debug_output_file"
} >&2
exit 1
    # Process files
    for file in "${fastq_files[@]}"; do
        local basename=$(basename "$file")
        local id=$(echo "$basename" | grep -oP "${PROJECT_CONFIG[BMC_FASTQ_ID_PATTERN]}\K[^_]*")
        local output_file="$output_dir/${id}${PROJECT_CONFIG[FASTQ_SUFFIX]}"
        
        log_info "Processing: $basename -> $(basename "$output_file")" "$log_file"
        
        if ! cat "$file" >> "$output_file"; then
            log_error "Failed to process file: $file" "$log_file"
            return 1
        fi
    done
    
    return 0
}
