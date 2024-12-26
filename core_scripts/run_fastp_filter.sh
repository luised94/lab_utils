#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=50G
#SBATCH --nice=10000
#SBATCH --exclude=c[5-22]
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luised94@mit.edu
# Script: run_fastp_filter.sh
# Purpose: Executes fastp filter script logic.
# Usage: sbatch --array=1-N%16 $HOME/lab_utils/core_scripts/run_fastp_filter.sh <experiment_directory>
# Author: Luis
# Date: 2024-12-26

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
readonly TOOL_DIR="${MONTH_DIR}/bowtie2_alignment"
readonly JOB_LOG_DIR="${TOOL_DIR}/job_${SLURM_ARRAY_JOB_ID}"
readonly TASK_LOG_DIR="${JOB_LOG_DIR}/task_${SLURM_ARRAY_TASK_ID}"
readonly TIMESTAMP=$(date +%Y%m%d_%H%M%S)
readonly MAIN_LOG="${TASK_LOG_DIR}/main_${TIMESTAMP}.log"
readonly ERROR_LOG="${TASK_LOG_DIR}/error_${TIMESTAMP}.log"
readonly PERFORMANCE_LOG="${TASK_LOG_DIR}/performance_${TIMESTAMP}.log"

# Create log directories
mkdir -p "${TASK_LOG_DIR}"
mkdir -p "${EXPERIMENT_DIR}/alignment"

# Log script start
log_message "INFO" "Starting alignment process for experiment: ${EXPERIMENT_DIR}"
log_message "INFO" "Job ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
log_message "INFO" "Log directory: ${TASK_LOG_DIR}"

# Load required modules
module purge
module load fastp/0.20.0

# Find fastq files
FASTQ_DIR="${EXPERIMENT_DIR}/fastq"
mapfile -t FASTQ_FILES < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name "consolidated*.fastq" | sort)
TOTAL_FILES=${#FASTQ_FILES[@]}

if [ $TOTAL_FILES -eq 0 ]; then
    log_message "ERROR" "No fastq files found in ${FASTQ_DIR}"
    exit 1
fi

# Get current fastq file
FASTQ_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
FASTQ_PATH="${FASTQ_FILES[$FASTQ_INDEX]}"

if [ -z "$FASTQ_PATH" ]; then
    log_message "Error: No fastq file found for index $FASTQ_INDEX"
    exit 1
fi

# Extract experiment ID from directory path
EXPERIMENT_ID=$(basename "$EXPERIMENT_DIR")

# Set quality parameters based on experiment ID
if [[ "$EXPERIMENT_ID" == "100303Bel" ]]; then
    # Legacy experiment from 2010
    QUALITY_THRESHOLD=20
    LENGTH_REQUIRED=20
    log_message "INFO" "Using legacy parameters for 2010 experiment: Q${QUALITY_THRESHOLD}, L${LENGTH_REQUIRED}"
else
    # Modern experiments
    QUALITY_THRESHOLD=30
    LENGTH_REQUIRED=50
    log_message "INFO" "Using standard parameters: Q${QUALITY_THRESHOLD}, L${LENGTH_REQUIRED}"
fi

# Extract and validate sample ID from filename
FASTQ_BASENAME=$(basename "$FASTQ_PATH")
if [[ ! "$FASTQ_BASENAME" =~ ^consolidate_.*_sequence\.fastq$ ]]; then
    log_message "ERROR" "Invalid input filename format: ${FASTQ_BASENAME}"
    log_message "ERROR" "Expected format: consolidate_<ID>_sequence.fastq"
    exit 1
fi

# Extract ID using sed with validation
SAMPLE_ID=$(echo "$FASTQ_BASENAME" | sed -n 's/consolidate_\(.*\)_sequence\.fastq/\1/p')
if [ -z "$SAMPLE_ID" ]; then
    log_message "ERROR" "Failed to extract sample ID from filename: ${FASTQ_BASENAME}"
    exit 1
fi

# Construct output path
OUTPUT_FASTQ="${EXPERIMENT_DIR}/fastq/processed_${SAMPLE_ID}_sequence.fastq"

# Validate output path doesn't already exist (optional)
if [ -f "$OUTPUT_FASTQ" ]; then
    log_message "WARNING" "Output file already exists: ${OUTPUT_FASTQ}"
    exit 1
fi

log_message "INFO" "Processing sample: ${SAMPLE_NAME}"
log_message "INFO" "Input: ${FASTQ_PATH}"
log_message "INFO" "Output: ${OUTPUT_FASTQ}"

# Alignment and sorting
log_message "INFO" "Starting fastp filtering"

if measure_performance "fastp_filtering" \
    fastp \
        --in1 "$FASTQ_PATH" \
        --out1 "$OUTPUT_FASTQ" \
        --cut_window_size 4 \
        --cut_mean_quality 20 \
        --n_base_limit 5 \
        --average_qual "$QUALITY_THRESHOLD" \
        --qualified_quality_phred "$QUALITY_THRESHOLD" \
        --unqualified_percent_limit 50 \
        --length_required "$LENGTH_REQUIRED" \
        --thread "$SLURM_CPUS_PER_TASK" \
        --compression 0 \
        --json "${TASK_LOG_DIR}/${SAMPLE_NAME}_fastp.json" \
        --html /dev/null; then
    
    log_message "INFO" "Performed fastp filtering."
else
    log_message "ERROR" "Fastp filtering failed for ${SAMPLE_NAME}"
    exit 1
fi
