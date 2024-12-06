#!/bin/bash
# functions moved

#SBATCH parameters from SLURM_CONFIG

set -euo pipefail

source "../functions/r_integration.sh"
source "../functions/bam_comparer.sh"

function main() {
    if [ $# -ne 2 ]; then
        log_error "Usage: $0 <directory> <time_id>"
        exit 1
    }
    
    local exp_dir="$HOME/data/$1"
    local time_id="$2"
    
    log_info "Starting BAM comparison"
    log_info "Job ID: ${SLURM_JOB_ID}"
    log_info "Array Task: ${SLURM_ARRAY_TASK_ID}"
    
    # Load required modules
    module purge
    for module in "${REQUIRED_MODULES[@]}"; do
        if ! module load "$module"; then
            log_error "Failed to load module: $module"
            exit 1
        fi
    done
    
    # Get BAM pairs from R
    mapfile -t bam_pairs < <(get_bam_pairs "$(basename "$exp_dir")" "$SLURM_ARRAY_TASK_ID")
    
    validate_bam_pairs "${bam_pairs[@]}" || exit 1
    
    local output_file=$(generate_output_name "${bam_pairs[0]}" \
                                           "${bam_pairs[1]}" \
                                           "$time_id" \
                                           "${exp_dir}/${OUTPUT_DIRS[BIGWIG]}")
    
    local threads=$((SLURM_CPUS_PER_TASK / 2))
    
    run_comparison "${bam_pairs[0]}" "${bam_pairs[1]}" "$output_file" "$threads" || exit 1
    
    log_info "Comparison completed successfully"
}

main "$@"
