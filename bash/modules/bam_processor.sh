#!/bin/bash

source "../config/quality_control_config.sh"

function setup_qc_directories() {
    local base_dir="$1"
    
    log_info "Setting up QC directories"
    
    for dir in "${QC_DIRS[@]}"; do
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
    
    log_info "Finding BAM files"
    
    local bam_files=$(find "$base_dir" -type f -name "*.bam")
    
    if [ -z "$bam_files" ]; then
        log_error "No BAM files found in: $base_dir"
        return 1
    }
    
    echo "$bam_files"
}

function get_output_basename() {
    local bam_path="$1"
    
    echo "$(basename "${bam_path%.bam}")"
}

function run_flagstat() {
    local bam_file="$1"
    local output_file="$2"
    
    log_info "Running flagstat on: $bam_file"
    
    if ! samtools flagstat ${SAMTOOLS_PARAMS[FLAGSTAT]} "$bam_file" > "$output_file"; then
        log_error "Flagstat failed for: $bam_file"
        return 1
    }
    
    return 0
}

function run_quickcheck() {
    local bam_file="$1"
    local output_file="$2"
    
    log_info "Running quickcheck on: $bam_file"
    
    if samtools quickcheck "$bam_file"; then
        echo -e 'QUICKCHECK\tTRUE' > "$output_file"
    else
        echo -e 'QUICKCHECK\tFALSE' > "$output_file"
        log_warning "Quickcheck failed for: $bam_file"
    fi
    
    return 0
}

function run_stats() {
    local bam_file="$1"
    local output_file="$2"
    
    log_info "Running stats on: $bam_file"
    
    if ! samtools stats ${SAMTOOLS_PARAMS[STATS]} "$bam_file" > "$output_file"; then
        log_error "Stats failed for: $bam_file"
        return 1
    }
    
    return 0
}

function process_bam_file() {
    local bam_file="$1"
    local qc_dir="$2"
    
    local basename=$(get_output_basename "$bam_file")
    local success=true
    
    # Run QC tools
    for tool in "${!BAM_QC_OUTPUTS[@]}"; do
        local output_file="${qc_dir}/${basename}${BAM_QC_OUTPUTS[$tool]}"
        
        case "$tool" in
            "FLAGSTAT")
                run_flagstat "$bam_file" "$output_file" || success=false
                ;;
            "QUICKCHECK")
                run_quickcheck "$bam_file" "$output_file" || success=false
                ;;
            "STATS")
                run_stats "$bam_file" "$output_file" || success=false
                ;;
        esac
    done
    
    $success
}
