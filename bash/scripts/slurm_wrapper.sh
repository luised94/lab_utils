#!/bin/bash
# functions moved

set -euo pipefail

source "../functions/slurm_job_handler.sh"
source "../functions/log_manager.sh"

function show_usage() {
    cat << EOF
Usage: $(basename "$0") <array_range> <script_name> <directory>

Arguments:
    array_range    SLURM array specification (e.g., "1-N%16")
    script_name    Script to execute
    directory      Experiment directory (without path)

Array Range Format:
    - Single task: "1"
    - Range: "1-10"
    - Range with limit: "1-100%16"
    - List: "1,2,5,10"

Example:
    $(basename "$0") "1-10%16" "align_fastq.sh" "240304Bel"
EOF
}

function main() {
    if [ $# -ne 3 ]; then
        show_usage
        exit 1
    }
    
    local array_range="$1"
    local script_name="$2"
    local dir_name="$3"
    
    validate_array_range "$array_range" || exit 1
    
    local script_path=$(find_script "$script_name") || exit 1
    local exp_dir=$(validate_experiment_dir "$dir_name") || exit 1
    
    local time_id=$(date "+${TIME_FORMATS[JOB_ID]}")
    
    log_info "Submitting SLURM job"
    log_info "Directory: $exp_dir"
    log_info "Script: $script_path"
    
    local job_id=$(sbatch --parsable \
                         --array="$array_range" \
                         "$script_path" \
                         "${exp_dir##*/}/" \
                         "$time_id")
    
    if [ -z "$job_id" ]; then
        log_error "Job submission failed"
        exit 1
    }
    
    log_info "Job submitted successfully"
    log_info "Job ID: $job_id"
    log_info "Time ID: $time_id"
    
    generate_log_commands "$exp_dir" "$job_id" "$time_id"
}

main "$@"
