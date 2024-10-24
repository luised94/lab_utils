#!/bin/bash
# functions moved

#SBATCH parameters from SLURM_CONFIG

set -euo pipefail

source "../functions/alignment_handler.sh"

function setup_directories() {
    local base_dir="$1"
    
    local dirs=(
        "logs"
        "alignment"
    )
    
    for dir in "${dirs[@]}"; do
        local full_path="${base_dir}/${dir}"
        mkdir -p "$full_path" || {
            log_error "Failed to create directory: $full_path"
            return 1
        }
    done
}

function main() {
    if [ $# -ne 2 ]; then
        log_error "Usage: $0 <experiment_dir> <time_id>"
        exit 1
    }
    
    local exp_dir="$HOME/data/$1"
    local time_id="$2"
    
    setup_directories "$exp_dir" || exit 1
    
    log_info "Starting alignment process"
    log_info "Job ID: ${SLURM_JOB_ID}"
    log_info "Array Task: ${SLURM_ARRAY_TASK_ID}"
    
    # Load required modules
    module purge
    for module in "${MODULE_REQUIREMENTS[@]}"; do
        if ! module load "$module"; then
            log_error "Failed to load module: $module"
            exit 1
        fi
    done
    
    # Get file lists
    mapfile -t fastq_files < <(find "$exp_dir" -type f -name "${FILE_PATTERNS[FASTQ]}")
    mapfile -t genome_files < <(find "$REFGENOME_DIR" -type f -name "${FILE_PATTERNS[GENOME]}")
    
    if [ ${#fastq_files[@]} -eq 0 ] || [ ${#genome_files[@]} -eq 0 ]; then
        log_error "No input files found"
        exit 1
    }
    
    log_info "Found ${#fastq_files[@]} FASTQ files and ${#genome_files[@]} genomes"
    
    # Calculate indices
    local indices=$(calculate_indices "$SLURM_ARRAY_TASK_ID" "${#fastq_files[@]}")
    local genome_index=${indices%:*}
    local fastq_index=${indices#*:}
    
    perform_alignment "${genome_files[$genome_index]}" \
                     "${fastq_files[$fastq_index]}" \
                     "${exp_dir}/alignment" \
                     "$SLURM_CPUS_PER_TASK" || exit 1
    
    log_info "Alignment task completed successfully"
}

main "$@"
