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
# Usage: sbatch --array=1-N%16 run_bamcoverage_array.sbatch <experiment_directory>
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

# Validate SLURM_ARRAY_TASK_ID
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Error: This script must be run as a SLURM array job"
    echo "Use: sbatch --array=1-N%16 $0 <experiment_directory>"
    exit 1
fi

# Logging functions
source "$HOME/lab_utils/core_scripts/functions_for_logging.sh"
readonly TOOL_NAME="blacklist_filtering"
eval "$(setup_logging ${TOOL_NAME})"

# Create log directories
mkdir -p "${EXPERIMENT_DIR}/alignment"

THREADS=$( nproc )

# Check if blacklist file exists
BLACKLIST_BED_FILE="$HOME/data/feature_files/20250423_merged_saccharomyces_cerevisiae_s288c_blacklist.bed"
if [[ ! -f "$BLACKLIST_BED_FILE" ]];
then
    echo "Blacklist file not found: $BLACKLIST_BED_FILE. Blacklist runs will be skipped."
    exit 1
fi

# Log script start
log_message "INFO" "Starting blacklist filtering process for experiment: ${EXPERIMENT_DIR}"
log_message "INFO" "Job ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
log_message "INFO" "Log directory: ${TASK_LOG_DIR}"

# Load required modules
module purge
module load python/2.7.13
module load deeptools/3.0.1
module load samtools
log_message "INFO" "Loaded required modules"

# Find BAM files
BAM_DIR="${EXPERIMENT_DIR}/alignment"
mapfile -t BAM_FILES < <(find "$BAM_DIR" -maxdepth 1 -type f -name "processed*_sorted.bam" | sort)
TOTAL_FILES=${#BAM_FILES[@]}
if [ "$TOTAL_FILES" -eq 0 ]; then
    log_message "ERROR" "No BAM files found in ${BAM_DIR}"
    exit 1
fi

# Validate array range
if [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL_FILES" ]; then
    log_message "WARNING" "Task ID ${SLURM_ARRAY_TASK_ID} exceeds number of jobs ${TOTAL_FILES}"
    exit 1
fi

# Get current BAM file
# Calculate array indices
BAM_INDEX=$((SLURM_ARRAY_TASK_ID - 1))

# Get current BAM file and normalization method
BAM_PATH="${BAM_FILES[$BAM_INDEX]}"
if [ ! -f "$BAM_PATH" ]; then
    log_message "WARNING" "Task ID ${SLURM_ARRAY_TASK_ID} bam path does not exist."
    log_message "WARNING" "Input: ${BAM_PATH}"
    exit 1
fi
# Set output name
SAMPLE_NAME=$(basename --suffix=_sorted.bam "$BAM_PATH" )
OUTPUT_BAM="${EXPERIMENT_DIR}/alignment/${SAMPLE_NAME}_blFiltered.bam"

log_message "INFO" "Processing sample: ${SAMPLE_NAME}"
log_message "INFO" "Input: ${BAM_PATH}"
log_message "INFO" "Output: ${OUTPUT_BAM}"

# Execute blacklist filtering ------
# Check if BAM file already exists
if [ -f "$OUTPUT_BAM" ]; then
  log_message "SKIP" "File already exists: ${OUTPUT_BAM}"

  if [ -f "${OUTPUT_BAM}.bai" ]; then
    log_message "SKIP" "File already exists: ${OUTPUT_BAM}.bai"
    exit 0
  fi

  if ! measure_performance "index" samtools index "$OUTPUT_BAM"; then
    log_message "ERROR" "BAM indexing failed for ${SAMPLE_NAME}"
    exit 1
  fi

  log_message "INFO" "Successfully completed processing for ${SAMPLE_NAME}"
  exit 0
fi

# Perform blacklist filtering
log_message "INFO" "Starting blacklist filtering processing"
if ! measure_performance "blacklist_filtering" \
  alignmentSieve --bam "${BAM_PATH}" \
                 --blackListFileName "$BLACKLIST_BED_FILE" \
                 --outFile "$OUTPUT_BAM" \
                 --numberOfProcessors $(( THREADS / 2 )); then
  log_message "ERROR" "Blacklist processing failed for ${SAMPLE_NAME}"
  exit 1
fi
log_message "INFO" "Successfully completed processing for ${SAMPLE_NAME}"

# Index the resulting BAM file
if ! measure_performance "index" samtools index "$OUTPUT_BAM"; then
  log_message "ERROR" "BAM indexing failed for ${SAMPLE_NAME}"
  exit 1
fi
log_message "INFO" "Successfully indexed for ${SAMPLE_NAME}"

# Log completion
log_message "INFO" "Task completed successfully"
