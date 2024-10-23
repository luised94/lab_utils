#!/bin/bash

source "../config/ngs_config.sh"

function validate_input() {
    local experiment_name="$1"
    
    log_info "Starting input validation"
    
    if [ $# -ne 1 ]; then
        log_error "Experiment name is required."
        return 1
    }
    
    if [[ "$experiment_name" == */ ]]; then
        log_error "Experiment name should not have trailing slash"
        return 1
    }
    
    local experiment_dir="${PATHS[BASE_DATA]}/$experiment_name"
    if [ ! -d "$experiment_dir" ]; then
        log_error "Experiment directory not found: $experiment_dir"
        return 1
    }
    
    echo "$experiment_dir"
}

function setup_output_directory() {
    local experiment_dir="$1"
    local output_dir="${experiment_dir}/${PATHS[FASTQ_SUBDIR]}"
    
    log_info "Creating output directory: $output_dir"
    mkdir -p "$output_dir"
    
    echo "$output_dir"
}

function find_fastq_files() {
    local experiment_dir="$1"
    local exclude_pattern="${NGS_CONFIG[EXCLUDE_PATTERNS]}"
    
    log_info "Searching for FASTQ files"
    
    find "$experiment_dir" -type f -name "*${PATTERNS[FASTQ_EXT]}" \
        ! \( -name "*unmapped*" -o -name "processed_*" \) | sort
}

function extract_unique_ids() {
    local -n files=$1
    local pattern="${PATTERNS[ID_PATTERN]}"
    
    log_info "Extracting unique IDs"
    printf '%s\n' "${files[@]}" | awk -F"$pattern" '{print $3}' | sort -u
}

function process_fastq_file() {
    local id="$1"
    local output_dir="$2"
    local -n files=$3
    
    local output_file="${output_dir}/${NGS_CONFIG[FILE_PREFIX]}${id}${NGS_CONFIG[FILE_SUFFIX]}"
    log_info "Processing ID: $id -> $output_file"
    
    local success=true
    for file in "${files[@]}"; do
        if [[ $file =~ $id ]]; then
            if ! cat "$file" >> "$output_file"; then
                log_error "Failed to append $file"
                success=false
                break
            fi
            if ! rm "$file"; then
                log_error "Failed to remove $file"
                success=false
                break
            fi
            log_info "Processed: $file"
        fi
    done
    
    $success
}

function verify_processing() {
    local experiment_dir="$1"
    local experiment_name="$2"
    local initial_count="$3"
    local processed_ids="$4"
    
    local current_count=$(find "$experiment_dir" -type f -name "*${PATTERNS[FASTQ_EXT]}" ! \( -name "*unmapped*" -o -name "processed_*" \) | wc -l)
    local sample_count=$(($(wc -l < "${experiment_dir}/${PATHS[DOC_SUBDIR]}/${experiment_name}_bmc_table.tsv") - 1))
    
    log_info "Processing Summary:"
    log_info "Initial files: $initial_count"
    log_info "Processed IDs: $processed_ids"
    log_info "Remaining files: $current_count"
    log_info "Expected samples: $sample_count"
    
    if [ "$processed_ids" -ne "$sample_count" ]; then
        log_warning "Mismatch between processed IDs ($processed_ids) and expected samples ($sample_count)"
    }
}
