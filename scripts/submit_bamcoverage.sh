#!/bin/bash
# submit_bamcoverage.sh

EXPERIMENT_DIR="$1"
if [ -z "$EXPERIMENT_DIR" ]; then
    echo "Error: Experiment directory not provided"
    echo "Usage: $0 <experiment_directory>"
    exit 1
fi

# Count BAM files
BAM_COUNT=$(find "${EXPERIMENT_DIR}/alignment" -maxdepth 1 -type f -name "*.sorted.bam" | wc -l)
echo "Found ${BAM_COUNT} BAM files"

if [ $BAM_COUNT -eq 0 ]; then
    echo "Error: No BAM files found in ${EXPERIMENT_DIR}/alignment"
    exit 1
fi

# Submit job
sbatch --array=1-${BAM_COUNT}%16 run_bamcoverage_array.sh "$EXPERIMENT_DIR"