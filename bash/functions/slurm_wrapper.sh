#!/bin/bash
# bash/functions/slurm_wrapper.sh

source "$HOME/lab_utils/bash/config/slurm_config.sh"
source "$HOME/lab_utils/bash/functions/logging_utils.sh"
source "$HOME/lab_utils/bash/functions/slurm_validator.sh"

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
#' Submit SLURM Array Job
#' @param array_range Character SLURM array specification
#' @param script_name Character Script name
#' @param dir_name Character Directory name
#' @param job_type Character Job type
#' @param log_file Character Log file path
#' @return Integer 0 if successful
submit_slurm_array_job() {
    local array_range="$1"
    local script_name="$2"
    local dir_name="$3"
    local job_type="$4"
    local timestamp_for_data="$5"
    local log_file="$6"
    
    local slurm_opts
    slurm_opts=$(get_slurm_options "$job_type")
    
    local job_id
    job_id=$(sbatch --parsable \
                    --array="$array_range" \
                    $slurm_opts \
                    "$script_name" \
                    "$dir_name"
                    "$timestamp")
    
    if [[ -z "$job_id" ]]; then
        log_error "Job submission failed" "$log_file"
        return 1
    fi
    
    echo "$job_id"
    return 0
}


#' Print Job Information
#' @param job_id Character SLURM job ID
#' @param log_file Character Log file path
print_job_info() {
    local job_id="$1"
    local log_file="$2"
    
    cat << EOF

Job Information:
  Job ID: $job_id
  Log File: $log_file

Monitor Commands:
  squeue -j $job_id
  sacct -j $job_id
  tail -f $log_file
  
Status Check:
  squeue -u $USER
  scontrol show job $job_id
EOF
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
    }
    
    # Format range
    local formatted_range
    if ! formatted_range=$(format_array_range "$range_input"); then
        log_error "Failed to format array range" "$log_file"
        return 1
    fi
    
    echo "$formatted_range"
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
    fi
    
    if [[ ! -x "$script_path" ]]; then
        log_error "Script not executable: $script_path" "$log_file"
        return 1
    fi
    
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
    fi
    
    echo "$exp_dir"
}