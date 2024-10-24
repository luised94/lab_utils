#!/bin/bash

#SBATCH parameters from SLURM_CONFIG

set -euo pipefail

source "../functions/fastq_processor.sh"

function main() {
    if [ $# -ne 2 ]; then
        log_error "Usage: $0 <directory> <time_id>"
        exit 1
    }
    
    local exp_dir="$HOME/data/$1"
    local time_id="$2"
    
    log_info "Starting FASTQ processing"
    log_info "Job ID: ${SLURM_JOB_ID}"
    log_info "Array Task: ${SLURM_ARRAY_TASK_ID}"
    
    setup_output_directories "$exp_dir" || exit 1
    
    # Load required modules
    module purge
    for module in "${REQUIRED_MODULES[@]}"; do
        if ! module load "$module"; then
            log_error "Failed to load module: $module"
            exit 1
        fi
    done
    
    # Get input files
    mapfile -t input_files < <(find_input_files "$exp_dir")
    
    if [ ${#input_files[@]} -eq 0 ]; then
        log_error "No input files found"
        exit 1
    }
    
    local task_index=$((SLURM_ARRAY_TASK_ID - 1))
    if [ $task_index -ge ${#input_files[@]} ]; then
        log_error "Task ID exceeds number of files"
        exit 1
    }
    
    process_fastq_file "${input_files[$task_index]}" "$exp_dir" || exit 1
    
    log_info "Processing completed successfully"
}

main "$@"
