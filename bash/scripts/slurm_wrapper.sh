#!/bin/bash
# bash/functions/slurm_wrapper.sh

source "$HOME/lab_utils/bash/config/slurm_config.sh"
source "$HOME/lab_utils/bash/functions/logging_utils.sh"

#' Validate SLURM Array Range
#' @param array_range Character SLURM array specification
#' @param log_file Character Log file path
#' @return Integer 0 if valid
validate_array_range() {
    local array_range="$1"
    local log_file="$2"
    
    # MIT-specific validation (must use %16)
    if [[ "$array_range" =~ ^[0-9]+-[0-9]+$ ]]; then
        log_error "Array range must include %16 limit: $array_range%16" "$log_file"
        return 1
    fi
    
    if [[ ! "$array_range" =~ ^([0-9]+|[0-9]+-[0-9]+%16|[0-9]+(,[0-9]+)*|[0-9]+-[0-9]+)$ ]]; then
        log_error "Invalid array range format: $array_range" "$log_file"
        return 1
    }
    
    return 0
}

#' Find Script in Project
#' @param script_name Character Script name
#' @param log_file Character Log file path
#' @return String Script path
find_slurm_script() {
    local script_name="$1"
    local log_file="$2"
    
    local script_path="$HOME/lab_utils/bash/scripts/$script_name"
    
    if [[ ! -f "$script_path" ]]; then
        log_error "Script not found: $script_path" "$log_file"
        return 1
    }
    
    if [[ ! -x "$script_path" ]]; then
        log_error "Script not executable: $script_path" "$log_file"
        return 1
    }
    
    echo "$script_path"
}

#' Validate Experiment Directory
#' @param dir_name Character Directory name
#' @param log_file Character Log file path
#' @return String Full directory path
validate_experiment_dir() {
    local dir_name="$1"
    local log_file="$2"
    
    local exp_dir="$HOME/data/$dir_name"
    
    if [[ ! -d "$exp_dir" ]]; then
        log_error "Experiment directory not found: $exp_dir" "$log_file"
        return 1
    }
    
    echo "$exp_dir"
}

#' Submit SLURM Job Main Function
#' @param array_range Character SLURM array specification
#' @param script_name Character Script name
#' @param dir_name Character Directory name
#' @return Integer 0 if successful
submit_slurm_job_main() {
    if [[ $# -ne 3 ]]; then
        show_usage
        return 1
    }
    
    # Initialize logging with timestamp for concurrent access
    local timestamp=$(date +%s)
    local log_file
    log_file=$(initialize_logging "slurm_submit_${timestamp}")
    
    local array_range="$1"
    local script_name="$2"
    local dir_name="$3"
    
    # Validate inputs
    if ! validate_array_range "$array_range" "$log_file"; then
        return 1
    fi
    
    local script_path
    if ! script_path=$(find_slurm_script "$script_name" "$log_file"); then
        return 1
    fi
    
    local exp_dir
    if ! exp_dir=$(validate_experiment_dir "$dir_name" "$log_file"); then
        return 1
    }
    
    # Submit job
    log_info "Submitting SLURM job:" "$log_file"
    log_info "  Directory: $exp_dir" "$log_file"
    log_info "  Script: $script_path" "$log_file"
    log_info "  Array: $array_range" "$log_file"
    
    local job_id
    job_id=$(sbatch --parsable --array="$array_range" "$script_path" "${dir_name}/" "$timestamp")
    
    if [[ -z "$job_id" ]]; then
        log_error "Job submission failed" "$log_file"
        return 1
    fi
    
    log_info "Job submitted successfully:" "$log_file"
    log_info "  Job ID: $job_id" "$log_file"
    log_info "  Time ID: $timestamp" "$log_file"
    
    # Output log locations
    cat << EOF

Job Information:
  Job ID: $job_id
  Time ID: $timestamp
  Log File: $log_file

Monitor Job:
  squeue -j $job_id
  sacct -j $job_id
  tail -f $log_file
EOF
    
    return 0
}

show_usage() {
    cat << EOF
Usage: $(basename "$0") <array_range> <script_name> <directory>

Arguments:
    array_range    SLURM array specification (must include %16)
    script_name    Script to execute
    directory      Experiment directory name

Examples:
    $(basename "$0") "1-10%16" "align_fastq.sh" "240304Bel"
    $(basename "$0") "1" "align_fastq.sh" "240304Bel"
    $(basename "$0") "1,2,5%16" "align_fastq.sh" "240304Bel"
EOF
}

# Execute if run as script
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    submit_slurm_job_main "$@"
fi
