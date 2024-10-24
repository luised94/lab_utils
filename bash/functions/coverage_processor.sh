#!/bin/bash

source "../config/coverage_config.sh"

function setup_coverage_directories() {
    local base_dir="$1"
    
    log_info "Setting up coverage directories"
    
    for dir in "${OUTPUT_DIRS[@]}"; do
        local full_path="$base_dir/$dir"
        mkdir -p "$full_path" || {
            log_error "Failed to create directory: $full_path"
            return 1
        }
    done
    
    return 0
}

function find_bam_files() {
    local base_dir="$1"
    local pattern="${FILE_PATTERNS[BAM_SUFFIX]}"
    
    log_info "Finding BAM files"
    
    local bam_files=$(find "$base_dir" -type f -name "*${pattern}" | sort)
    
    if [ -z "$bam_files" ]; then
        log_error "No BAM files found matching pattern: $pattern"
        return 1
    }
    
    echo "$bam_files"
}

function generate_output_name() {
    local bam_path="$1"
    local time_id="$2"
    local output_dir="$3"
    
    local basename=$(basename "${bam_path%.bam}")
    echo "${output_dir}/${time_id}_${basename}${FILE_PATTERNS[BIGWIG_SUFFIX]}"
}

function build_bamcoverage_command() {
    local input_file="$1"
    local output_file="$2"
    
    local cmd="bamCoverage"
    cmd+=" -b \"$input_file\""
    cmd+=" -o \"$output_file\""
    cmd+=" --binSize ${COVERAGE_PARAMS[BIN_SIZE]}"
    cmd+=" --normalizeUsing ${COVERAGE_PARAMS[NORMALIZATION]}"
    cmd+=" --minMappingQuality ${COVERAGE_PARAMS[MIN_MAPPING_QUALITY]}"
    
    if [ "${COVERAGE_PARAMS[IGNORE_DUPLICATES]}" = true ]; then
        cmd+=" --ignoreDuplicates"
    fi
    
    echo "$cmd"
}

function process_bam_coverage() {
    local bam_file="$1"
    local output_file="$2"
    
    log_info "Processing coverage for: $bam_file"
    log_info "Output: $output_file"
    
    local cmd=$(build_bamcoverage_command "$bam_file" "$output_file")
    
    log_info "Executing: $cmd"
    if ! eval "$cmd"; then
        log_error "Coverage generation failed for: $bam_file"
        return 1
    }
    
    log_info "Successfully generated coverage for: $bam_file"
    return 0
}
