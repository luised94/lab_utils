#!/bin/bash

source "../config/bam_comparison_config.sh"

function get_bam_pairs() {
    local base_dir="$1"
    local task_id="$2"
    local r_script="${R_CONFIG[SCRIPT_PATH]}"
    
    log_info "Getting BAM pairs from R script"
    
    if [ ! -f "$r_script" ]; then
        log_error "R script not found: $r_script"
        return 1
    fi
    
    local pairs=$(Rscript "$r_script" "$base_dir" "$task_id" | grep -v '^NULL$')
    
    if [ -z "$pairs" ]; then
        log_error "No output from R script"
        return 1
    fi
    
    echo "$pairs"
}

function validate_bam_pairs() {
    local -a pairs=("$@")
    
    log_info "Validating BAM pairs"
    
    if [ ${#pairs[@]} -ne 2 ]; then
        log_error "Expected 2 BAM files, got ${#pairs[@]}"
        return 1
    fi
    
    for bam in "${pairs[@]}"; do
        if [ ! -f "$bam" ]; then
            log_error "BAM file not found: $bam"
            return 1
        fi
    done
    
    return 0
}
