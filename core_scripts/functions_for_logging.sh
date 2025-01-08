#!/bin/bash
# Purpose: Contains simple reusable functions for the logging portion of slurm scripts.
# Usage: source ~/lab_utils/core_scripts/functions_for_logging.sh
# Author: Luis
# Date: 2025-01-08

# Function to display usage
display_usage() {
    echo "Usage: sbatch --array=1-N%16 $0 <experiment_directory>"
    echo "Example: sbatch --array=1-10%16 $0 /home/user/data/240304Bel"
    echo "Note: Array range should not exceed the number of fastq files"
    exit 1
}

# Function to log messages
log_message() {
    local level=$1
    local message=$2
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] [${level}] [Task ${SLURM_ARRAY_TASK_ID}] ${message}" | tee -a "${MAIN_LOG}"
}

# Function to log performance metrics
log_performance() {
    local stage=$1
    local duration=$2
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] ${stage}: ${duration} seconds" >> "${PERFORMANCE_LOG}"
}

# Function to measure command execution time
measure_performance() {
    local stage=$1
    shift
    local start_time=$(date +%s)
    "$@" 2>> "${ERROR_LOG}"
    local status=$?
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    log_performance "${stage}" "${duration}"
    return $status
}

# Function to validate array range
validate_array_range() {
    local total_files=$1
    local array_start=$(echo $SLURM_ARRAY_TASK_MIN)
    local array_end=$(echo $SLURM_ARRAY_TASK_MAX)
    log_message "Validating array range..."
    log_message "Total fastq files: $total_files"
    log_message "Array range: $array_start-$array_end"
    if [ $array_end -gt $total_files ]; then
        log_message "WARNING: Array range ($array_end) exceeds number of fastq files ($total_files)"
        log_message "Suggestion: Use --array=1-${total_files}%16"
    fi
}

setup_logging() {
    local tool_name="$1"

    if [[ -z "$tool_name" || ! "$tool_name" =~ ^[a-zA-Z0-9_-]+$ ]]; then
        if [[ -z "$tool_name" ]]; then
            echo "Error: Tool name cannot be empty" >&2
        else
            echo "Error: Invalid tool name. Must contain only alphanumeric characters, underscores, and hyphens." >&2
        fi
        return 1
    fi

    if [[ ${#tool_name} -gt 50 ]]; then
        echo "Error: Tool name too long (max 50 characters)" >&2
        return 1
    fi

    local log_root="$HOME/logs" # Or get from config/env var
    local current_month=$(date +%Y-%m)
    local month_dir="${log_root}/${current_month}"
    local tool_dir="${month_dir}/${tool_name}"
    local job_log_dir="${tool_dir}/job_${SLURM_ARRAY_JOB_ID}"
    local task_log_dir="${job_log_dir}/task_${SLURM_ARRAY_TASK_ID}"
    local timestamp=$(date +%Y%m%d_%H%M%S)

    if ! mkdir -p "$task_log_dir"; then
        echo "Error: Failed to create log directory: $task_log_dir" >&2
        return 1
    fi

    printf "MAIN_LOG='%s'\nERROR_LOG='%s'\nPERFORMANCE_LOG='%s'\nTASK_LOG_DIR='%s'\nJOB_LOG_DIR='%s'\nTOOL_DIR='%s'\n" \
           "${task_log_dir}/main_${timestamp}.log" "${task_log_dir}/error_${timestamp}.log" "${task_log_dir}/performance_${timestamp}.log" \
           "${task_log_dir}" "${job_log_dir}" "${tool_dir}"
}
