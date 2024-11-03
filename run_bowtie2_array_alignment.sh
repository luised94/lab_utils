#!/bin/bash

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
EXPERIMENT_DIR="$1"

# Validate SLURM_ARRAY_TASK_ID
if [ -z "$SLURM_ARRAY_TASK_ID" ]; then
    echo "Error: This script must be run as a SLURM array job"
    echo "Use: sbatch --array=1-N%16 $0 <experiment_directory>"
    exit 1
fi

# Constants
GENOME_DIR="$HOME/data/REFGENS/SaccharomycescerevisiaeS288C"
GENOME_INDEX="$GENOME_DIR/SaccharomycescerevisiaeS288C_index"
LOG_DIR="$HOME/logs"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_FILE="${LOG_DIR}/ngs_alignment_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${TIMESTAMP}.log"

# Create required directories
mkdir -p "${LOG_DIR}"
mkdir -p "${EXPERIMENT_DIR}/alignment"

# Function to log messages
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
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

# Function to measure command execution time
measure_performance() {
    local start_time=$(date +%s)
    "$@"
    local status=$?
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    log_message "Command: $1"
    log_message "Duration: ${duration} seconds"
    log_message "Exit Status: ${status}"
    return $status
}

# Validate directories and files
if [ ! -d "$EXPERIMENT_DIR" ]; then
    log_message "Error: Experiment directory not found: $EXPERIMENT_DIR"
    exit 1
fi

if [ ! -f "${GENOME_INDEX}.1.bt2" ]; then
    log_message "Error: Genome index not found: $GENOME_INDEX"
    exit 1
fi

# Load required modules
module purge
module load bowtie2
module load samtools

log_message "Starting alignment process"
log_message "Experiment Directory: $EXPERIMENT_DIR"
log_message "Array Task ID: $SLURM_ARRAY_TASK_ID"

# Find fastq files
FASTQ_DIR="${EXPERIMENT_DIR}/fastq"
mapfile -t FASTQ_FILES < <(find "$FASTQ_DIR" -maxdepth 1 -type f -name "*.fastq")
TOTAL_FILES=${#FASTQ_FILES[@]}

if [ $TOTAL_FILES -eq 0 ]; then
    log_message "Error: No fastq files found in $FASTQ_DIR"
    exit 1
fi

# Validate array range against number of files
validate_array_range $TOTAL_FILES

# Get current fastq file
FASTQ_INDEX=$((SLURM_ARRAY_TASK_ID - 1))
FASTQ_PATH="${FASTQ_FILES[$FASTQ_INDEX]}"

if [ -z "$FASTQ_PATH" ]; then
    log_message "Error: No fastq file found for index $FASTQ_INDEX"
    exit 1
fi

# Generate output name
SAMPLE_NAME=$(basename "$FASTQ_PATH" .fastq)
OUTPUT_BAM="${EXPERIMENT_DIR}/alignment/${SAMPLE_NAME}.sorted.bam"

log_message "Processing sample: $SAMPLE_NAME"
log_message "Input FASTQ: $FASTQ_PATH"
log_message "Output BAM: $OUTPUT_BAM"

# Alignment and sorting
log_message "Starting alignment and sorting"
if measure_performance bowtie2 -x "$GENOME_INDEX" \
                              -U "$FASTQ_PATH" \
                              -p "$SLURM_CPUS_PER_TASK" 2>> "$LOG_FILE" | \
   samtools view -@ "$SLURM_CPUS_PER_TASK" -bS - 2>> "$LOG_FILE" | \
   samtools sort -@ "$SLURM_CPUS_PER_TASK" -o "$OUTPUT_BAM" - 2>> "$LOG_FILE"; then
    
    # Index BAM file
    log_message "Starting BAM indexing"
    if measure_performance samtools index "$OUTPUT_BAM" 2>> "$LOG_FILE"; then
        log_message "Successfully completed processing for $SAMPLE_NAME"
    else
        log_message "Error: BAM indexing failed for $SAMPLE_NAME"
        exit 1
    fi
else
    log_message "Error: Alignment/sorting failed for $SAMPLE_NAME"
    exit 1
fi
