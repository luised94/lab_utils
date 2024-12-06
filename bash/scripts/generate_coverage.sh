#!/bin/bash
# functions moved

#SBATCH parameters from SLURM_CONFIG

set -euo pipefail

source "../functions/coverage_processor.sh"

function main() {
    if [ $# -ne 2 ]; then
        log_error "Usage: $0 <directory> <time_id>"
        exit 1
    }
    
    local exp_dir="$HOME/data/$1"
    local time_id="$2"
    
    log_info "Starting coverage generation"
    log_info "Job ID: ${SLURM_JOB_ID}"
    log_info "Array Task: ${SLURM_ARRAY_TASK_ID}"
    
    setup_coverage_directories "$exp_dir" || exit 1
    
    # Load required modules
    module purge
    for module in "${REQUIRED_MODULES[@]}"; do
        if ! module load "$module"; then
            log_error "Failed to load module: $module"
            exit 1
        fi
    done
    
    # Get BAM files
    mapfile -t bam_files < <(find_bam_files "$exp_dir")
    
    if [ ${#bam_files[@]} -eq 0 ]; then
        log_error "No BAM files found"
        exit 1
    }
    
    local task_index=$((SLURM_ARRAY_TASK_ID - 1))
    if [ $task_index -ge ${#bam_files[@]} ]; then
        log_error "Task ID exceeds number of files"
        exit 1
    }
    
    local output_file=$(generate_output_name "${bam_files[$task_index]}" \
                                           "$time_id" \
                                           "${exp_dir}/${OUTPUT_DIRS[BIGWIG]}")
    
    process_bam_coverage "${bam_files[$task_index]}" "$output_file" || exit 1
    
    log_info "Coverage generation completed successfully"
}

main "$@"
