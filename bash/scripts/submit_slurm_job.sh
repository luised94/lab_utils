#!/bin/bash

source "$HOME/lab_utils/bash/functions/slurm_wrapper.sh"
#' Submit SLURM Job Main Function
#' @param array_range Character SLURM array specification
#' @param script_name Character Script name
#' @param dir_name Character Directory name
#' @return Integer 0 if successful
submit_slurm_job_main() {
    if [[ $# -ne 3 ]]; then
        show_usage
        return 1
    fi
    
    # Initialize logging with timestamp for concurrent access
    local timestamp_for_log=$(date +%s)
    local log_file
    log_file=$(initialize_logging "slurm_submit_${timestamp_for_log}")
    
    local array_range="$1"
    local script_name="$2"
    local dir_name="$3"
    local timestamp_for_data=$(date "+%Y-%m-%d-%H-%M-%S")
    
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
    fi
    
    # Submit job
    log_info "Submitting SLURM job:" "$log_file"
    log_info "  Directory: $exp_dir" "$log_file"
    log_info "  Script: $script_path" "$log_file"
    log_info "  Array: $array_range" "$log_file"
    
    local job_id
    if ! job_id=$(submit_slurm_array_job "$array_range" "$script_name" "$exp_dir" "$job_type"  "$timestamp_for_data" "$log_file"); then
        return 1
    fi
    
    if [[ -z "$job_id" ]]; then
        log_error "Job submission failed" "$log_file"
        return 1
    fi
    
    print_job_info "$job_id"  "$log_file"
    return 0
}


show_usage() {
    cat << EOF
Usage: $(basename "$0") <array_range> <script_name> <directory>

Arguments:
    array_range    SLURM array specification:
                   - Single task: "1"
                   - Multiple tasks: "1,2,5"
                   - Range: "1-10" (automatically adds %16)
    script_name    Script to execute
    directory      Experiment directory name

Examples:
    $(basename "$0") "1-10" "align_fastq.sh" "240304Bel"
    $(basename "$0") "1" "align_fastq.sh" "240304Bel"
    $(basename "$0") "1,2,5" "align_fastq.sh" "240304Bel"

Note: submit_slurm_job automatically applies %16 limit to ranges
EOF
}

# Execute if run as script
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    submit_slurm_job_main "$@"
fi
