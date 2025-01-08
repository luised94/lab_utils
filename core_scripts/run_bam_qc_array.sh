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















# Quality Control Configuration
readonly MIN_MAPPING_QUALITY=20
readonly MIN_INSERT_SIZE=0
readonly MAX_INSERT_SIZE=1000

# Create log and output directories
mkdir -p "${TASK_LOG_DIR}"
readonly QUALITY_CONTROL_DIR="${EXPERIMENT_DIR}/quality_control/bam"
mkdir -p "${QUALITY_CONTROL_DIR}"

# Logging functions
source $HOME/lab_utils/core_scripts/functions_for_logging.sh
readonly TOOL_NAME="REPLACE_ME"
eval "$(setup_logging ${TOOL_NAME})"







# Performance logging functions







# Measure command execution time












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
