#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=50G
#SBATCH --nice=10000
#SBATCH --exclude=c[5-22]
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luised94@mit.edu
# Script: run_macs2_peak_calling.sh
# Purpose: Executes macs2 peak calling via the miniforge installation.
# Usage: sbatch --array=1-N%16 run_bamcoverage_array.sh <experiment_directory>
# Dependencies: submit_macs2_peak_calling.sh, macs2, miniforge, slurm
# Date: 2024-12-26

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

# Logging setup
readonly CURRENT_MONTH=$(date +%Y-%m)
readonly LOG_ROOT="$HOME/logs"
readonly MONTH_DIR="${LOG_ROOT}/${CURRENT_MONTH}"
readonly TOOL_DIR="${MONTH_DIR}/peak_calling"
readonly JOB_LOG_DIR="${TOOL_DIR}/job_${SLURM_ARRAY_JOB_ID}"
readonly TASK_LOG_DIR="${JOB_LOG_DIR}/task_${SLURM_ARRAY_TASK_ID}"
readonly TIMESTAMP=$(date +%Y%m%d_%H%M%S)
readonly MAIN_LOG="${TASK_LOG_DIR}/main_${TIMESTAMP}.log"
readonly ERROR_LOG="${TASK_LOG_DIR}/error_${TIMESTAMP}.log"
readonly PERFORMANCE_LOG="${TASK_LOG_DIR}/performance_${TIMESTAMP}.log"

# Create log directories
mkdir -p "${TASK_LOG_DIR}"
mkdir -p "${EXPERIMENT_DIR}/peak_calling"


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
# Source setup
if ! source $HOME/lab_utils/core_scripts/setup_conda_and_macs2.sh; then
    log_message "ERROR" "Failed to activate MACS2 environment"
    exit 1
fi

# Print environment info for debugging
log_message "INFO" "Active conda env: $CONDA_DEFAULT_ENV"
macs_path=$(which macs2)
log_message "INFO" "MACS2 path: $macs_path"

# Find BAM files
BAM_DIR="${EXPERIMENT_DIR}/alignment"
mapfile -t BAM_FILES < <(find "$BAM_DIR" -maxdepth 1 -type f -name "*_sorted.bam" | sort)
TOTAL_FILES=${#BAM_FILES[@]}
if [ $TOTAL_FILES -eq 0 ]; then
    log_message "ERROR" "No BAM files found in ${BAM_DIR}"
    exit 1
fi

# Validate array range
if [ $SLURM_ARRAY_TASK_ID -gt $TOTAL_FILES ]; then
    log_message "WARNING" "Task ID ${SLURM_ARRAY_TASK_ID} exceeds number of jobs ${TOTAL_FILES}"
    exit 1
fi

# Get current BAM file
# Calculate array indices
BAM_INDEX=$(( SLURM_ARRAY_TASK_ID - 1))

# Get current BAM file and normalization method
BAM_PATH="${BAM_FILES[$BAM_INDEX]}"
if [ ! -f $BAM_PATH ]; then
    log_message "WARNING" "Task ID ${SLURM_ARRAY_TASK_ID} bam path does not exist."
    log_message "WARNING" "Input: ${BAM_PATH}"
    exit 1
fi

# Set output name
SAMPLE_NAME=$(basename --suffix=_sorted.bam "$BAM_PATH" )
OUTPUT_PEAK="${EXPERIMENT_DIR}/peak_calling/${SAMPLE_NAME}_"

log_message "INFO" "Processing sample: ${SAMPLE_NAME}"
log_message "INFO" "Normalization method: ${NORM_METHOD}"
log_message "INFO" "Input: ${BAM_PATH}"
log_message "INFO" "Output: ${OUTPUT_PEAK}"

# Parameter explanation:
# --nomodel : Skip the model building step
# --extsize 147 : Typical yeast nucleosome size (~150bp)
# --shift -73 : Half the fragment size (147/2)
# --mfold 3 100 : Minimum fold enrichment of 3
# --pvalue 1e-6 : Matches paper's significance threshold
# --bdg : Generate bedGraph files for visualization
# --SPMR : Normalize to signal per million reads

# Execute bamCoverage
log_message "INFO" "Starting bamCoverage processing"
if measure_performance "peak_calling" \
    macs2 callpeak \
        -t "$TEST_SAMPLE" \
        -c "$INPUT_CONTROL" \
        -n "${OUTPUT_PREFIX}_241010Bel" \
        -g "1.2e7" \
        --nomodel \
        --extsize 147 \
        --shift -73 \
        --keep-dup auto \
        --pvalue 1e-6 \
        --mfold 3 100 \
        --outdir "$OUTDIR" \
        --bdg \
        --SPMR; then
    log_message "INFO" "Successfully completed processing for ${SAMPLE_NAME}"
else
    log_error "bamCoverage processing failed for ${SAMPLE_NAME}"
    exit 1
fi

# Log completion
log_message "INFO" "Task completed successfully"
