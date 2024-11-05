#!/bin/bash
# run_fastqc_array.sh
# Purpose: Execute FastQC as a SLURM array job
# Version: 1.1.0
# Compatibility: Bash 4.2+, SLURM

# Strict error handling
set -euo pipefail

# Validate input arguments
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <experiment_directory> <fastq_subdirectory>"
    exit 1
fi

# Parse arguments
readonly EXPERIMENT_DIR="$1"
readonly FASTQ_SUBDIR="$2"

# Validate SLURM array job context
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "Error: This script must be run as a SLURM array job"
    exit 1
fi

# Logging configuration
readonly CURRENT_MONTH=$(date +%Y-%m)
readonly LOG_ROOT="${HOME}/logs"
readonly MONTH_DIR="${LOG_ROOT}/${CURRENT_MONTH}"
readonly TOOL_DIR="${MONTH_DIR}/fastqc"
readonly JOB_LOG_DIR="${TOOL_DIR}/job_${SLURM_ARRAY_JOB_ID}"
readonly TASK_LOG_DIR="${JOB_LOG_DIR}/task_${SLURM_ARRAY_TASK_ID}"
readonly TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Log file paths
readonly MAIN_LOG="${TASK_LOG_DIR}/main_${TIMESTAMP}.log"
readonly ERROR_LOG="${TASK_LOG_DIR}/error_${TIMESTAMP}.log"
readonly PERFORMANCE_LOG="${TASK_LOG_DIR}/performance_${TIMESTAMP}.log"

# Create log directories
mkdir -p "${TASK_LOG_DIR}"
readonly QUALITY_CONTROL_DIR="${EXPERIMENT_DIR}/quality_control"
mkdir -p "${QUALITY_CONTROL_DIR}"

# Logging functions
log_message() {
    local level="$1"
    local message="$2"
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

# Find FASTQ files
FASTQ_DIR="${EXPERIMENT_DIR}/${FASTQ_SUBDIR}"
mapfile -t FASTQ_FILES < <(find "$FASTQ_DIR" -maxdepth 1 \( -name "*.fastq" -o -name "*.fq" \))
TOTAL_FILES=${#FASTQ_FILES[@]}

# Validate array task
if [[ ${SLURM_ARRAY_TASK_ID} -gt ${TOTAL_FILES} ]]; then
    log_message "ERROR" "Task ID exceeds number of files"
    exit 1
fi

# Select current file
FASTQ_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
FASTQ_PATH="${FASTQ_FILES[${FASTQ_INDEX}]}"
SAMPLE_NAME=$(basename "${FASTQ_PATH}" | sed -E 's/\.(fastq|fq)$//')

# Load required modules
module purge
module load fastqc
module load java

# Execute FastQC with performance measurement
log_message "INFO" "Processing sample: ${SAMPLE_NAME}"
log_message "INFO" "Input file: ${FASTQ_PATH}"

if measure_performance "fastqc" fastqc \
    --outdir "${QUALITY_CONTROL_DIR}" \
    --threads "${SLURM_CPUS_PER_TASK:-1}" \
    "${FASTQ_PATH}"; then
    
    log_message "INFO" "FastQC completed successfully for ${SAMPLE_NAME}"
else
    log_message "ERROR" "FastQC processing failed for ${SAMPLE_NAME}"
    exit 1
fi

# Log completion
log_message "INFO" "Task completed successfully"
