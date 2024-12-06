#!/bin/bash
# run_bam_qc_array.sh
# Purpose: Execute BAM Quality Control using Samtools
# Version: 1.1.0
# Compatibility: Bash 4.2+, SLURM

# Strict error handling
set -euo pipefail

# Validate input arguments
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <experiment_directory> <bam_subdirectory>"
    exit 1
fi

# Parse arguments
readonly EXPERIMENT_DIR="$1"
readonly BAM_SUBDIR="$2"

# Validate SLURM array job context
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "Error: This script must be run as a SLURM array job"
    exit 1
fi

# Logging configuration
readonly CURRENT_MONTH=$(date +%Y-%m)
readonly LOG_ROOT="${HOME}/logs"
readonly MONTH_DIR="${LOG_ROOT}/${CURRENT_MONTH}"
readonly TOOL_DIR="${MONTH_DIR}/bam_qc"
readonly JOB_LOG_DIR="${TOOL_DIR}/job_${SLURM_ARRAY_JOB_ID}"
readonly TASK_LOG_DIR="${JOB_LOG_DIR}/task_${SLURM_ARRAY_TASK_ID}"
readonly TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Log file paths
readonly MAIN_LOG="${TASK_LOG_DIR}/main_${TIMESTAMP}.log"
readonly ERROR_LOG="${TASK_LOG_DIR}/error_${TIMESTAMP}.log"
readonly PERFORMANCE_LOG="${TASK_LOG_DIR}/performance_${TIMESTAMP}.log"

# Quality Control Configuration
readonly MIN_MAPPING_QUALITY=20
readonly MIN_INSERT_SIZE=0
readonly MAX_INSERT_SIZE=1000

# Create log and output directories
mkdir -p "${TASK_LOG_DIR}"
readonly QUALITY_CONTROL_DIR="${EXPERIMENT_DIR}/quality_control/bam"
mkdir -p "${QUALITY_CONTROL_DIR}"

# Logging functions
log_message() {
    local level="$1"
    local message="$2"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] [${level}] [Task ${SLURM_ARRAY_TASK_ID}] ${message}" | tee -a "${MAIN_LOG}"
}

# Performance logging functions
log_performance() {
    local stage="$1"
    local duration="$2"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] ${stage}: ${duration} seconds" >> "${PERFORMANCE_LOG}"
}

# Measure command execution time
measure_performance() {
    local stage="$1"
    shift
    local start_time=$(date +%s)
    "$@" 2>> "${ERROR_LOG}"
    local status=$?
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    log_performance "${stage}" "${duration}"
    return $status
}

# Find BAM files
BAM_DIR="${EXPERIMENT_DIR}/${BAM_SUBDIR}"
mapfile -t BAM_FILES < <(find "$BAM_DIR" -maxdepth 1 \( -name "*.sorted.bam" -o -name "*.bam" \))
TOTAL_FILES=${#BAM_FILES[@]}

# Validate array task
if [[ ${SLURM_ARRAY_TASK_ID} -gt ${TOTAL_FILES} ]]; then
    log_message "ERROR" "Task ID exceeds number of files"
    exit 1
fi

# Select current file
BAM_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
BAM_PATH="${BAM_FILES[${BAM_INDEX}]}"
SAMPLE_NAME=$(basename "${BAM_PATH}" | sed -E 's/\.(sorted\.)?bam$//')

# Load required modules
module purge
module load samtools

# Comprehensive Samtools Quality Control
log_message "INFO" "Generating comprehensive BAM quality control for ${SAMPLE_NAME}"

# 1. Basic Statistics
log_message "INFO" "Generating Samtools stats"
measure_performance "samtools_stats" samtools stats \
    -q "${MIN_MAPPING_QUALITY}" \
    "${BAM_PATH}" > "${QUALITY_CONTROL_DIR}/${SAMPLE_NAME}.samtools_stats.txt"

# 2. Flagstat (Quick mapping statistics)
log_message "INFO" "Generating Flagstat report"
measure_performance "samtools_flagstat" samtools flagstat \
    "${BAM_PATH}" > "${QUALITY_CONTROL_DIR}/${SAMPLE_NAME}.flagstat.txt"

# 3. Idxstats (Chromosome-level statistics)
log_message "INFO" "Generating Idxstats report"
measure_performance "samtools_idxstats" samtools idxstats \
    "${BAM_PATH}" > "${QUALITY_CONTROL_DIR}/${SAMPLE_NAME}.idxstats.txt"

# 4. Insert Size Distribution
log_message "INFO" "Calculating insert size distribution"
measure_performance "samtools_insert_size" samtools stats \
    -i "${MIN_INSERT_SIZE}:${MAX_INSERT_SIZE}" \
    "${BAM_PATH}" > "${QUALITY_CONTROL_DIR}/${SAMPLE_NAME}.insert_size.txt"

# 5. Coverage Depth
log_message "INFO" "Calculating coverage depth"
measure_performance "samtools_depth" samtools depth \
    -a "${BAM_PATH}" > "${QUALITY_CONTROL_DIR}/${SAMPLE_NAME}.depth.txt"

# 6. Generate a comprehensive report
log_message "INFO" "Generating comprehensive BAM report"
{
    echo "Sample: ${SAMPLE_NAME}"
    echo "BAM File: ${BAM_PATH}"
    echo "----------------------------"
    echo "Samtools Stats Summary:"
    grep "^SN" "${QUALITY_CONTROL_DIR}/${SAMPLE_NAME}.samtools_stats.txt" | \
        sed 's/^SN\t//'
} > "${QUALITY_CONTROL_DIR}/${SAMPLE_NAME}.bam_qc_summary.txt"

# Log completion
log_message "INFO" "BAM Quality Control completed for ${SAMPLE_NAME}"
