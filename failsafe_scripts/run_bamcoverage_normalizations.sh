#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=50G
#SBATCH --nice=10000
#SBATCH --exclude=c[5-22]
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luised94@mit.edu
# Script: run_bamcoverage_array.sh
# Purpose: Executes deepTools bamCoverage as SLURM array job for multiple BAM files
# Usage: sbatch --array=1-N%16 run_bamcoverage_array.sh <experiment_directory>
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

# Normalization methods array
declare -a NORM_METHODS=("RPKM" "CPM" "BPM" "RPGC")

# coverage parameters - hardcoded values
BIN_SIZE=10
EFFECTIVE_GENOME_SIZE=12157105
MIN_MAPPING_QUALITY=20

# Logging setup
readonly CURRENT_MONTH=$(date +%Y-%m)
readonly LOG_ROOT="$HOME/logs"
readonly MONTH_DIR="${LOG_ROOT}/${CURRENT_MONTH}"
readonly TOOL_DIR="${MONTH_DIR}/bamcoverage"
readonly JOB_LOG_DIR="${TOOL_DIR}/job_${SLURM_ARRAY_JOB_ID}"
readonly TASK_LOG_DIR="${JOB_LOG_DIR}/task_${SLURM_ARRAY_TASK_ID}"
readonly TIMESTAMP=$(date +%Y%m%d_%H%M%S)
readonly MAIN_LOG="${TASK_LOG_DIR}/main_${TIMESTAMP}.log"
readonly ERROR_LOG="${TASK_LOG_DIR}/error_${TIMESTAMP}.log"
readonly PERFORMANCE_LOG="${TASK_LOG_DIR}/performance_${TIMESTAMP}.log"

# Create log directories
mkdir -p "${TASK_LOG_DIR}"
mkdir -p "${EXPERIMENT_DIR}/coverage"


# Logging functions
log_message() {
    local level=$1
    local message=$2
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] [${level}] [Task ${SLURM_ARRAY_TASK_ID}] ${message}" | tee -a "${MAIN_LOG}"
}

log_error() {
    log_message "ERROR" "$1" >&2
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" >> "${ERROR_LOG}"
}

log_performance() {
    local stage=$1
    local duration=$2
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo "[${timestamp}] ${stage}: ${duration} seconds" >> "${PERFORMANCE_LOG}"
}

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

# Log script start
log_message "INFO" "Starting bamCoverage process for experiment: ${EXPERIMENT_DIR}"
log_message "INFO" "Job ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
log_message "INFO" "Log directory: ${TASK_LOG_DIR}"

# Load required modules
module purge
module load python
module load deeptools

# Find BAM files
BAM_DIR="${EXPERIMENT_DIR}/alignment"
mapfile -t BAM_FILES < <(find "$BAM_DIR" -maxdepth 1 -type f -name "*_sorted.bam" | sort)
TOTAL_FILES=${#BAM_FILES[@]}
TOTAL_JOBS=$((TOTAL_FILES * ${#NORM_METHODS[@]}))
if [ $TOTAL_FILES -eq 0 ]; then
    log_message "ERROR" "No BAM files found in ${BAM_DIR}"
    exit 1
fi

# Validate array range
if [ $SLURM_ARRAY_TASK_ID -gt $TOTAL_JOBS ]; then
    log_message "WARNING" "Task ID ${SLURM_ARRAY_TASK_ID} exceeds number of jobs ${TOTAL_FILES}"
    exit 1
fi

# Get current BAM file
# Calculate array indices
BAM_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) / ${#NORM_METHODS[@]} ))
NORM_INDEX=$(( (SLURM_ARRAY_TASK_ID - 1) % ${#NORM_METHODS[@]} ))

# Get current BAM file and normalization method
BAM_PATH="${BAM_FILES[$BAM_INDEX]}"
NORM_METHOD="${NORM_METHODS[$NORM_INDEX]}"
if [ ! -f $BAM_PATH ]; then
    log_message "WARNING" "Task ID ${SLURM_ARRAY_TASK_ID} bam path does not exist."
    log_message "WARNING" "Input: ${BAM_PATH}"
    exit 1
fi
# Set output name
SAMPLE_NAME=$(basename "$BAM_PATH" --suffix="_sorted.bam")
OUTPUT_BIGWIG="${EXPERIMENT_DIR}/coverage/${SAMPLE_NAME}_${NORM_METHOD}.bw"

log_message "INFO" "Processing sample: ${SAMPLE_NAME}"
log_message "INFO" "Normalization method: ${NORM_METHOD}"
log_message "INFO" "Input: ${BAM_PATH}"
log_message "INFO" "Output: ${OUTPUT_BIGWIG}"

# Build bamCoverage command with specific parameters for each method
COMMON_PARAMS="--bam ${BAM_PATH} \
    --outFileName ${OUTPUT_BIGWIG} \
    --binSize ${BIN_SIZE} \
    --minMappingQuality ${MIN_MAPPING_QUALITY} \
    --ignoreDuplicates \
    --normalizeUsing ${NORM_METHOD} \
    --numberOfProcessors ${SLURM_CPUS_PER_TASK}"

# Add RPGC-specific parameters
if [ "${NORM_METHOD}" == "RPGC" ]; then
    COMMON_PARAMS+=" --effectiveGenomeSize ${EFFECTIVE_GENOME_SIZE}"
    log_message "INFO" "Added effectiveGenomeSize parameter for RPGC normalization"
fi

# Execute bamCoverage
log_message "INFO" "Starting bamCoverage processing"
if measure_performance "bamcoverage" bamCoverage $COMMON_PARAMS; then
    log_message "INFO" "Successfully completed processing for ${SAMPLE_NAME}"
else
    log_error "bamCoverage processing failed for ${SAMPLE_NAME}"
    exit 1
fi

log_message "INFO" "Parameters used:\n"
log_message "INFO" "$COMMON_PARAMS"
# Log completion
log_message "INFO" "Task completed successfully"
