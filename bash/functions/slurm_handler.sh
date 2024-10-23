#!/bin/bash

source "../config/genome_index_config.sh"

function setup_logging() {
    local log_dir="${GENOME_CONFIG[LOG_DIR]}"
    local timestamp=$(date "+%Y-%m-%d-%H-%M-%S")
    local job_id="${SLURM_ARRAY_JOB_ID:-standalone}"
    
    mkdir -p "$log_dir" || {
        echo "Failed to create log directory: $log_dir"
        return 1
    }
    
    echo "${log_dir}/indexing_${job_id}_${timestamp}"
}

function validate_slurm_env() {
    if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
        log_error "This script must be run as a SLURM array job"
        return 1
    }
    
    if [ -z "${SLURM_ARRAY_JOB_ID:-}" ]; then
        log_error "SLURM_ARRAY_JOB_ID not set"
        return 1
    }
    
    return 0
}

function load_required_modules() {
    log_info "Loading required modules"
    
    module purge || log_warning "Failed to purge modules"
    
    for module in "${MODULES[@]}"; do
        if ! module load "$module"; then
            log_error "Failed to load module: $module"
            return 1
        fi
        log_info "Loaded module: $module"
    done
    
    return 0
}
