#!/bin/bash
# submit_alignment.sh

EXPERIMENT_DIR="$1"
if [ -z "$EXPERIMENT_DIR" ]; then
    echo "Error: Experiment directory not provided"
    echo "Usage: $0 <experiment_directory>"
    exit 1
fi

# Count fastq files
FASTQ_COUNT=$(find "${EXPERIMENT_DIR}/fastq" -maxdepth 1 -type f -name "*.fastq" | wc -l)
echo "Found ${FASTQ_COUNT} fastq files"

if [ $FASTQ_COUNT -eq 0 ]; then
    echo "Error: No fastq files found in ${EXPERIMENT_DIR}/fastq"
    exit 1
fi

# Submit job
sbatch --array=1-${FASTQ_COUNT}%16 run_bowtie2_array_alignment.sh "$EXPERIMENT_DIR"
