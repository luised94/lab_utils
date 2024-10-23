#!/bin/bash

source "../config/slurm_config.sh"

function validate_array_range() {
    local range="$1"
    
    log_info "Validating array range: $range"
    
    if [[ ! "$range" =~ ^[0-9]+(|-[0-9]+(%[0-9]+)?|,[0-9]+)*$ ]]; then
        log_error "Invalid array range format: $range"
        return 1
    }
    
    if [[ "$range" =~ %([0-9]+) ]]; then
        local concurrent="${BASH_REMATCH[1]}"
        if ((concurrent > ${SLURM_WRAPPER[MAX_ARRAY_SIZE]})); then
            log_warning "Concurrent jobs ($concurrent) exceeds recommended maximum (${SLURM_WRAPPER[MAX_ARRAY_SIZE]})"
        }
    }
    
    return 0
}

function find_script() {
    local script_name="$1"
    local base_dir="${SLURM_WRAPPER[SCRIPT_BASE_DIR]}"
    
    log_info "Searching for script: $script_name"
    
    local script_path=$(find "$base_dir" -type f -name "$script_name")
    
    if [ -z "$script_path" ]; then
        log_error "Script not found: $script_name"
        return 1
    }
    
    if [ ! -x "$script_path" ]; then
        log_error "Script not executable: $script_path"
        return 1
    }
    
    echo "$script_path"
}

function validate_experiment_dir() {
    local dir_name="$1"
    local base_dir="${SLURM_WRAPPER[DATA_BASE_DIR]}"
    
    log_info "Validating experiment directory: $dir_name"
    
    local full_path=$(find -H "$base_dir" -maxdepth 1 -type d -name "$dir_name")
    
    if [ -z "$full_path" ]; then
        log_error "Directory not found: $dir_name"
        return 1
    }
    
    echo "$full_path"
}
