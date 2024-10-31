#!/bin/bash
# bash/functions/slurm_wrapper.sh

source "$HOME/lab_utils/bash/config/slurm_config.sh"
source "$HOME/lab_utils/bash/functions/logging_utils.sh"
source "$HOME/lab_utils/bash/functions/slurm_validator.sh"

get_slurm_options() {
    local job_type="$1"
    # [existing implementation]
}

submit_slurm_array_job() {
    local array_range="$1"
    local script_name="$2"
    local dir_name="$3"
    local job_type="$4"
    local log_file="$5"
    # [existing implementation]
}

print_job_info() {
    local job_id="$1"
    local log_file="$2"
    # [existing implementation]
}

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
    fi
    
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
