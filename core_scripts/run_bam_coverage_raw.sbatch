#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=50G
#SBATCH --nice=10000
#SBATCH --exclude=c[5-22]
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luised94@mit.edu
# Script: run_bamcoverage_array.sbatch
# Purpose: Executes deepTools bamCoverage as SLURM array job for multiple BAM files
# Usage: sbatch --array=1-N%16 run_bam_coverage_raw.sbatch <experiment_directory>
# Date: 2024-11-03

# Function to display usage
display_usage() {
    echo "Usage: sbatch --array=1-N%16 $0 <experiment_directory>"
    echo "Example: sbatch --array=1-10%16 $0 /home/user/data/240304Bel"
    exit 1
}

# Validate input arguments
if [ "$#" -ne 1 ]; then
    display_usage
fi

# Parse arguments
EXPERIMENT_DIR="$1"
EXPERIMENT_ID=$( basename "$EXPERIMENT_DIR" )

# Validate SLURM_ARRAY_TASK_ID
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Error: This script must be run as a SLURM array job"
    echo "Use: sbatch --array=1-N%16 $0 <experiment_directory>"
    exit 1
fi

# Logging functions
source "$HOME/lab_utils/core_scripts/functions_for_logging.sh"
readonly TOOL_NAME="coverage"
eval "$(setup_logging "${EXPERIMENT_ID}" ${TOOL_NAME})"

# coverage parameters - hardcoded values
BIN_SIZE=10
MIN_MAPPING_QUALITY=20

# Create log directories
mkdir -p "${EXPERIMENT_DIR}/coverage"

# Log script start
log_message "INFO" "Starting bamCoverage process for experiment: ${EXPERIMENT_DIR}"
log_message "INFO" "Job ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
log_message "INFO" "Log directory: ${TASK_LOG_DIR}"

# Load required modules
module purge
module load python/2.7.13
module load deeptools/3.0.1

# Find BAM files
BAM_DIR="${EXPERIMENT_DIR}/alignment"
mapfile -t BAM_FILES < <(find "$BAM_DIR" -maxdepth 1 -type f -name "*.bam" | sort)
TOTAL_FILES=${#BAM_FILES[@]}
if [ "$TOTAL_FILES" -eq 0 ]; then
    log_message "ERROR" "No BAM files found in ${BAM_DIR}"
    exit 1
fi

# Validate array range
if [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL_JOBS" ]; then
    log_message "WARNING" "Task ID ${SLURM_ARRAY_TASK_ID} exceeds number of jobs ${TOTAL_FILES}"
    exit 1
fi

# Get current BAM file
# Calculate array indices
BAM_INDEX=$(( SLURM_ARRAY_TASK_ID - 1 ))

# Get current BAM file and normalization method
BAM_PATH="${BAM_FILES[$BAM_INDEX]}"
if [ ! -f "$BAM_PATH" ]; then
    log_message "WARNING" "Task ID ${SLURM_ARRAY_TASK_ID} bam path does not exist."
    log_message "WARNING" "Input: ${BAM_PATH}"
    exit 1
fi
# Set output name
SAMPLE_NAME=$(basename --suffix=.bam "$BAM_PATH" )
OUTPUT_BIGWIG="${EXPERIMENT_DIR}/coverage/${SAMPLE_NAME}_RAW.bw"

log_message "INFO" "Processing sample: ${SAMPLE_NAME}"
log_message "INFO" "Normalization method: RAW"
log_message "INFO" "Input: ${BAM_PATH}"
log_message "INFO" "Output: ${OUTPUT_BIGWIG}"

# Build bamCoverage command with specific parameters for each method
COMMON_PARAMS="--bam ${BAM_PATH} \
    --outFileName ${OUTPUT_BIGWIG} \
    --binSize ${BIN_SIZE} \
    --minMappingQuality ${MIN_MAPPING_QUALITY} \
    --numberOfProcessors ${SLURM_CPUS_PER_TASK}"

log_message "INFO" "Parameters used:\n"
log_message "INFO" "$COMMON_PARAMS"

if [ -f "$OUTPUT_BIGWIG" ]; then
  log_message "SKIP" "File already exists: ${OUTPUT_BIGWIG}"
  log_message "INFO" "Successfully completed processing for ${SAMPLE_NAME}"
  exit 0
fi

# Execute bamCoverage
log_message "INFO" "Starting bamCoverage processing"
# shellcheck disable=SC2086 # Intentional: $COMMON_PARAMS should be word-split
if ! measure_performance "bamcoverage" bamCoverage $COMMON_PARAMS; then
    log_error "bamCoverage processing failed for ${SAMPLE_NAME}"
    exit 1
fi

log_message "INFO" "Successfully completed processing for ${SAMPLE_NAME}"
# Log completion
log_message "INFO" "Task completed successfully"
