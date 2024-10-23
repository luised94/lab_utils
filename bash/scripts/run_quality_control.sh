#!/bin/bash

#SBATCH parameters from SLURM_CONFIG

set -euo pipefail

source "../functions/quality_control.sh"

function validate_input() {
    local dir_name="$1"
    
    if [ -z "$dir_name" ]; then
        log_error "Directory name not provided"
        return 1
    }
    
    local full_path="$HOME/data/$dir_name"
    if [ ! -d "$full_path" ]; then
        log_error "Directory not found: $full_path"
        return 1
    }
    
    echo "$full_path"
}

function main() {
    if [ $# -lt 1 ]; then
        log_error "Usage: $0 <directory_name>"
        exit 1
    }
    
    local base_dir=$(validate_input "$1") || exit 1
    local qc_dir=$(setup_qc_directory "$base_dir") || exit 1
    
    load_required_modules || exit 1
    
    mapfile -t raw_files < <(find_fastq_files "$base_dir")
    mapfile -t processed_files < <(find_processed_fastq "$base_dir")
    
    local task_index=$((SLURM_ARRAY_TASK_ID - 1))
    
    if [ -n "${raw_files[$task_index]:-}" ]; then
        run_fastqc "${raw_files[$task_index]}" "$qc_dir" || exit 1
    fi
    
    if [ -n "${processed_files[$task_index]:-}" ]; then
        run_fastqc "${processed_files[$task_index]}" "$qc_dir" || exit 1
    fi
    
    log_info "Quality control completed successfully"
}

main "$@"
