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
    fi
    
    # Check SLURM environment
    if [[ -z "${SLURM_JOB_ID:-}" ]]; then
        log_error "Not running in SLURM environment" "$log_file"
        return 1
    fi
    
    # Verify reference directory
    if [[ ! -d "${SLURM_GENOMES[BASE_DIR]}" ]]; then
        log_error "Reference genome directory not found: ${SLURM_GENOMES[BASE_DIR]}" "$log_file"
        return 1
    fi
    
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

#' Format Array Range for MIT SLURM
#' @param range_input String Input range (e.g., "1-10", "1,2,3", "5")
#' @return String Formatted range with %16
format_array_range() {
    local range_input="$1"
    
    # Single number
    if [[ "$range_input" =~ ^[0-9]+$ ]]; then
        echo "$range_input"
        return 0
    fi
    
    # Comma-separated list
    if [[ "$range_input" =~ ^[0-9]+(,[0-9]+)*$ ]]; then
        echo "$range_input"
        return 0
    fi
    
    # Range (add %16 if not present)
    if [[ "$range_input" =~ ^[0-9]+-[0-9]+$ ]]; then
        echo "${range_input}%16"
        return 0
    fi
    
    # Invalid format
    return 1
}
#' Validate SLURM Array Range
#' @param range_input String Input range
#' @param log_file Character Log file path
#' @return Integer 0 if valid
validate_array_range() {
    local range_input="$1"
    local log_file="$2"
    
    # Validate format
    if ! [[ "$range_input" =~ ^([0-9]+|[0-9]+-[0-9]+|[0-9]+(,[0-9]+)*)$ ]]; then
        log_error "Invalid array range format: $range_input" "$log_file"
        return 1
    fi
    
    # Format range
    local formatted_range
    if ! formatted_range=$(format_array_range "$range_input"); then
        log_error "Failed to format array range" "$log_file"
        return 1
    fi
    
    echo "$formatted_range"
    return 0
}


