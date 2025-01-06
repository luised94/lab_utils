#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=50G
#SBATCH --nice=10000
#SBATCH --exclude=c[5-22]
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luised94@mit.edu
# Script: run_bowtie2_array_alignment.sh
# Purpose: Executes bowtie2 alignment as SLURM array job for multiple fastq files
# Usage: sbatch --array=1-N%16 run_bowtie2_array_alignment.sh <experiment_directory>
# Author: [Your Name]
# Date: 2024-11-03

# Function to display usage
display_usage() {
    echo "Usage: sbatch --array=1-N%16 $0 <experiment_directory>"
    echo "Example: sbatch --array=1-10%16 $0 /home/user/data/240304Bel"
    echo "Note: Array range should not exceed the number of fastq files"
    exit 1
}

# Validate input arguments
if [ "$#" -ne 1 ]; then
    display_usage
fi

# Parse arguments
EXPERIMENT_DIR=$(realpath "$1")


# Validate SLURM_ARRAY_TASK_ID
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Error: This script must be run as a SLURM array job"
    echo "Use: sbatch --array=1-N%16 $0 <experiment_directory>"
    exit 1
fi

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

# Constants
GENOME_DIR="$HOME/data/REFGENS/SaccharomycescerevisiaeS288C"
GENOME_INDEX="$GENOME_DIR/SaccharomycescerevisiaeS288C_index"

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


if [ ! -f "${GENOME_INDEX}.1.bt2" ]; then
    log_message "Error: Genome index not found: $GENOME_INDEX"
    exit 1
fi

# Log script start
log_message "INFO" "Starting alignment process for experiment: ${EXPERIMENT_DIR}"
log_message "INFO" "Job ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
log_message "INFO" "Log directory: ${TASK_LOG_DIR}"

# Load required modules
module purge
module load bowtie2
module load samtools

# Find fastq files
FASTQ_DIR="${EXPERIMENT_DIR}/fastq"
mapfile -t FASTQ_FILES < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name "consolidated*.fastq" | sort)
#mapfile -t FASTQ_FILES < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name "processed*.fastq" | sort)
TOTAL_FILES=${#FASTQ_FILES[@]}

if [ $TOTAL_FILES -eq 0 ]; then
    log_message "ERROR" "No fastq files found in ${FASTQ_DIR}"
    exit 1
fi

# Validate array range against number of files
#validate_array_range $TOTAL_FILES

# Get current fastq file
FASTQ_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
FASTQ_PATH="${FASTQ_FILES[$FASTQ_INDEX]}"

if [ -z "$FASTQ_PATH" ]; then
    log_message "Error: No fastq file found for index $FASTQ_INDEX"
    exit 1
fi

# Generate output name
SAMPLE_NAME=$(basename --suffix=.fastq "$FASTQ_PATH" )
OUTPUT_BAM="${EXPERIMENT_DIR}/alignment/${SAMPLE_NAME}_to_S288C_sorted.bam"

log_message "INFO" "Processing sample: ${SAMPLE_NAME}"
log_message "INFO" "Input: ${FASTQ_PATH}"
log_message "INFO" "Output: ${OUTPUT_BAM}"

# Alignment and sorting
log_message "INFO" "Starting alignment and sorting"
if measure_performance "alignment_and_sorting" \
    bowtie2 -x "$GENOME_INDEX" \
            -U "$FASTQ_PATH" \
            -p "$SLURM_CPUS_PER_TASK" 2>> "${ERROR_LOG}" | \
    samtools view -@ "$SLURM_CPUS_PER_TASK" -bS - 2>> "${ERROR_LOG}" | \
    samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$OUTPUT_BAM" - 2>> "${ERROR_LOG}"; then
# "-q --mp 4 --met-stderr"

    log_message "INFO" "Starting BAM indexing"
    if measure_performance "indexing" samtools index "$OUTPUT_BAM"; then
        log_message "INFO" "Successfully completed processing for ${SAMPLE_NAME}"
    else
        log_message "ERROR" "BAM indexing failed for ${SAMPLE_NAME}"
        exit 1
    fi
else
    log_message "ERROR" "Alignment/sorting failed for ${SAMPLE_NAME}"
    exit 1
fi

# Log completion
log_message "INFO" "Task completed successfully"
