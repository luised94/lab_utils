#!/bin/bash
# functions moved

#SBATCH parameters from SLURM_CONFIG

set -euo pipefail

source "../functions/bam_processor.sh"

function main() {
    if [ $# -ne 2 ]; then
        log_error "Usage: $0 <directory> <time_id>"
        exit 1
    }
    
    local exp_dir="$HOME/data/$1"
    local time_id="$2"
    
    log_info "Starting BAM quality control"
    log_info "Job ID: ${SLURM_JOB_ID}"
    log_info "Array Task: ${SLURM_ARRAY_TASK_ID}"
    
    setup_qc_directories "$exp_dir" || exit 1
    
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
    
    process_bam_file "${bam_files[$task_index]}" "${exp_dir}/${QC_DIRS[OUTPUT]}" || exit 1
    
    log_info "Quality control completed successfully"
}

main "$@"