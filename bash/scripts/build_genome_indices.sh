#!/bin/bash

#SBATCH -N ${SLURM_CONFIG[NODES]}
#SBATCH -n ${SLURM_CONFIG[TASKS]}
#SBATCH --mem-per-cpu=${SLURM_CONFIG[MEM_PER_CPU]}
#SBATCH --exclude=${SLURM_CONFIG[EXCLUDE_NODES]}
#SBATCH --mail-type=${SLURM_CONFIG[MAIL_TYPE]}
#SBATCH --mail-user=${SLURM_CONFIG[MAIL_USER]}
# functions moved

set -euo pipefail

source "../functions/slurm_handler.sh"
source "../functions/genome_indexer.sh"

function main() {
    validate_slurm_env || exit 1
    
    local log_base=$(setup_logging) || exit 1
    exec 1>"${log_base}.out" 2>"${log_base}.err"
    
    log_info "Starting genome indexing job"
    log_info "Job ID: ${SLURM_JOB_ID}, Array ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
    
    load_required_modules || exit 1
    
    mapfile -t genome_paths < <(find_reference_genomes) || exit 1
    
    if [ ${#genome_paths[@]} -eq 0 ]; then
        log_error "No reference genomes found"
        exit 1
    }
    
    local task_index=$((SLURM_ARRAY_TASK_ID - 1))
    if [ $task_index -ge ${#genome_paths[@]} ]; then
        log_error "Task ID exceeds number of genomes"
        exit 1
    }
    
    build_genome_index "${genome_paths[$task_index]}" || exit 1
    
    log_info "Job completed successfully"
}

main
