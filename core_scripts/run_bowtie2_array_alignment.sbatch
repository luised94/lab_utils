#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=90G
#SBATCH --nice=10000
#SBATCH --exclude=c[5-22]
#SBATCH --mail-type=ALL
#SBATCH --mail-user=luised94@mit.edu
# Script: run_bowtie2_array_alignment.sbatch
# Purpose: Executes bowtie2 alignment as SLURM array job for multiple fastq files
# Usage: sbatch --array=1-N%16 run_bowtie2_array_alignment.sbatch <experiment_directory>
# Author: [Your Name]
# Date: 2024-11-03

# Validate input arguments
if [ "$#" -ne 1 ]; then
    echo "Usage: sbatch --array=1-N%16 $0 <experiment_directory>"
    echo "Example: sbatch --array=1-10%16 $0 /home/user/data/240304Bel"
    echo "Note: Array range should not exceed the number of fastq files"
    exit 1
fi

# Parse arguments
EXPERIMENT_DIR=$(realpath "$1")
EXPERIMENT_ID=$( basename "$EXPERIMENT_DIR" )

# Validate SLURM_ARRAY_TASK_ID
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Error: This script must be run as a SLURM array job"
    echo "Use: sbatch --array=1-N%16 $0 <experiment_directory>"
    exit 1
fi

# Function to log messages
source $HOME/lab_utils/core_scripts/functions_for_logging.sh
readonly TOOL_NAME="bowtie2"
eval "$(setup_logging "${EXPERIMENT_ID}" ${TOOL_NAME} ${SLURM_ARRAY_JOB_ID} ${SLURM_ARRAY_TASK_ID} )"

# Constants
GENOME_DIR="$HOME/data/REFGENS/SaccharomycescerevisiaeS288C"
GENOME_INDEX="$GENOME_DIR/SaccharomycescerevisiaeS288C_index"

mkdir -p "${EXPERIMENT_DIR}/alignment"

if [ ! -f "${GENOME_INDEX}.1.bt2" ]; then
    log_message "Error: Genome index not found: $GENOME_INDEX"
    exit 1
fi

# Log script start
log_message "INFO" "Starting alignment process for experiment: ${EXPERIMENT_DIR}"
log_message "INFO" "Job ID: ${SLURM_ARRAY_JOB_ID}, Task ID: ${SLURM_ARRAY_TASK_ID}"
log_message "INFO" "Log directory: ${LOG_DIR}"

# Load required modules
module purge
module load bowtie2
module load samtools

# Find fastq files
FASTQ_DIR="${EXPERIMENT_DIR}/fastq"
mapfile -t FASTQ_FILES < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name "processed*.fastq" | sort)
TOTAL_FILES=${#FASTQ_FILES[@]}

if [ $TOTAL_FILES -eq 0 ]; then
    log_message "ERROR" "No fastq files found in ${FASTQ_DIR}"
    exit 1
fi

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

if [ -f "${OUTPUT_BAM}" ]; then
  log_message "SKIP" "File already exists: ${OUTPUT_BAM}"
  log_message "INFO" "Successfully completed processing for ${SAMPLE_NAME}"
  exit 0
fi

log_message "INFO" "Processing sample: ${SAMPLE_NAME}"
log_message "INFO" "Input: ${FASTQ_PATH}"
log_message "INFO" "Output: ${OUTPUT_BAM}"

# Alignment and sorting
# Most sensitive alignment mode
# Report only the best alignment
# Limit seed attempts
# Match bonus
# Maximum penalty for mismatch
# Penalty for ambiguous bases
# Read gap open, extend penalties
# Reference gap open, extend penalties
# Alignment score threshold
#--rdg 5,1 \
#--rfg 5,1 \
#--score-min L,0,-0.3 \
#--ma 2 \
log_message "INFO" "Starting alignment and sorting"
if measure_performance "alignment_and_sorting" \
    bowtie2 -x "$GENOME_INDEX" \
            -U "$FASTQ_PATH" \
            -p "$SLURM_CPUS_PER_TASK" \
            --sensitive \
            -k 1 \
            -D 15 \
            -R 2 \
            --mp 4 \
            --np 1 \
            2>> "${ERROR_LOG}" | \
    samtools view -@ "$SLURM_CPUS_PER_TASK" -bS - 2>> "${ERROR_LOG}" | \
    samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$OUTPUT_BAM" - 2>> "${ERROR_LOG}"; then
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
