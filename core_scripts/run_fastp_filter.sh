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

# Validate input arguments
if [ "$#" -ne 1 ]; then
    display_usage
fi

# Function to log messages
source $HOME/lab_utils/core_scripts/functions_for_logging.sh
readonly TOOL_NAME="fastp"
eval "$(setup_logging ${TOOL_NAME})"

# Parse arguments
EXPERIMENT_DIR="$1"

# Validate SLURM_ARRAY_TASK_ID
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Error: This script must be run as a SLURM array job"
    echo "Use: sbatch --array=1-N%16 $0 <experiment_directory>"
    exit 1
fi

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

# Extract and validate sample ID from filename
FASTQ_BASENAME=$(basename "$FASTQ_PATH")
if [[ ! "$FASTQ_BASENAME" =~ ^consolidated_.*_sequence\.fastq$ ]]; then
    log_message "ERROR" "Invalid input filename format: ${FASTQ_BASENAME}"
    log_message "ERROR" "Expected format: consolidate_<ID>_sequence.fastq"
    exit 1
fi

# Extract ID using sed with validation
SAMPLE_ID=$(echo "$FASTQ_BASENAME" | sed -n 's/consolidated_\(.*\)_sequence\.fastq/\1/p')
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

log_message "INFO" "Processing sample: ${SAMPLE_ID}"
log_message "INFO" "Input: ${FASTQ_PATH}"
log_message "INFO" "Output: ${OUTPUT_FASTQ}"

# Alignment and sorting
log_message "INFO" "Starting fastp filtering"

# Extract experiment ID from directory path
EXPERIMENT_ID=$(basename "$EXPERIMENT_DIR")

# Set quality parameters based on experiment ID
if [[ "$EXPERIMENT_ID" == "100303Bel" ]]; then
    # Legacy experiment from 2010
    MINIMUM_BASE_QUALITY=20
    MINIMUM_READ_LENGTH=20
    log_message "INFO" "Using legacy parameters for 2010 experiment: Q${MINIMUM_BASE_QUALITY}, L${MINIMUM_READ_LENGTH}"
else
    # Modern experiments
    MINIMUM_BASE_QUALITY=20
    MINIMUM_READ_LENGTH=50
    log_message "INFO" "Using standard parameters: Q${MINIMUM_BASE_QUALITY}, L${MINIMUM_READ_LENGTH}"
fi

# Quality Filtering Parameters
readonly MAXIMUM_UNQUALIFIED_BASE_PERCENT=50
readonly MAXIMUM_N_BASE_COUNT=5

# Read Length Parameters
readonly MAXIMUM_READ_LENGTH=150

# Complexity Filtering
readonly COMPLEXITY_WINDOW_SIZE=4
#readonly COMPLEXITY_THRESHOLD=30
readonly OVERREPRESENTATION_SAMPLING=50
readonly DUPLICATION_CALC_ACCURACY=3

# Performance Parameters
readonly COMPRESSION_LEVEL=0
readonly CPU_THREADS="$SLURM_CPUS_PER_TASK"

# Options are not available in fastp 0.20.0 version available in the linux cluster.
#--dedup \
#--dup_calc_accuracy "$DUPLICATION_CALC_ACCURACY" \
if measure_performance "fastp_filtering" \
    fastp \
        --in1 "$FASTQ_PATH" \
        --out1 "$OUTPUT_FASTQ" \
        --adapter_sequence auto \
        --cut_window_size "$COMPLEXITY_WINDOW_SIZE" \
        --cut_mean_quality "$MINIMUM_BASE_QUALITY" \
        --cut_front \
        --cut_tail \
        --cut_right \
        --n_base_limit "$MAXIMUM_N_BASE_COUNT" \
        --average_qual "$MINIMUM_BASE_QUALITY" \
        --qualified_quality_phred "$MINIMUM_BASE_QUALITY" \
        --unqualified_percent_limit "$MAXIMUM_UNQUALIFIED_BASE_PERCENT" \
        --length_required "$MINIMUM_READ_LENGTH" \
        --thread "$CPU_THREADS" \
        --overrepresentation_analysis \
        --overrepresentation_sampling "$OVERREPRESENTATION_SAMPLING" \
        --compression "$COMPRESSION_LEVEL" \
        --json "${TASK_LOG_DIR}/${SAMPLE_ID}_fastp.json" \
        --html "${TASK_LOG_DIR}/${SAMPLE_ID}_fastp.html"; then
    
    log_message "INFO" "Performed fastp filtering."
else
    log_message "ERROR" "Fastp filtering failed for ${SAMPLE_ID}"
    exit 1
fi
