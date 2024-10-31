#!/bin/bash
# bash/functions/slurm_validator.sh

source "$HOME/lab_utils/bash/config/slurm_config.sh"
source "$HOME/lab_utils/bash/functions/logging_utils.sh"

#' Validate SLURM Environment
#' @param job_type Character Type of job (ALIGN, QC, BW)
#' @param log_file Character Log file path
#' @return Integer 0 if valid
validate_slurm_environment() {
    local job_type="$1"
    local log_file="$2"
    
    # Validate resource requirements
    if [[ ! " ALIGN QC BW " =~ " $job_type " ]]; then
        log_error "Invalid job type: $job_type" "$log_file"
        return 1
    }
    
    # Check SLURM environment
    if [[ -z "${SLURM_JOB_ID:-}" ]]; then
        log_error "Not running in SLURM environment" "$log_file"
        return 1
    }
    
    # Verify reference directory
    if [[ ! -d "${SLURM_GENOMES[BASE_DIR]}" ]]; then
        log_error "Reference genome directory not found: ${SLURM_GENOMES[BASE_DIR]}" "$log_file"
        return 1
    }
    
    return 0
}

#' Validate Module Availability
#' @param modules Array Required module names
#' @param log_file Character Log file path
#' @return Integer 0 if all modules available
validate_modules() {
    local -a modules=("$@")
    local log_file="${modules[-1]}" # Last argument is log file
    unset 'modules[-1]'            # Remove log file from array
    
    for module in "${modules[@]}"; do
        if ! module avail "$module" 2>/dev/null; then
            log_error "Required module not available: $module" "$log_file"
            return 1
        fi
    done
    
    return 0
}

validate_array_range() {
    local array_range="$1"
    local log_file="$2"
    # MIT-specific validation
    # [implementation from previous wrapper]
}
