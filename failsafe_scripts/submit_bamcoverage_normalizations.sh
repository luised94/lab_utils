#!/bin/bash
# submit_bamcoverage_normalizations.sh
EXPERIMENT_DIR="$1"
if [ -z "$EXPERIMENT_DIR" ]; then
    echo "Error: Experiment directory not provided"
    echo "Usage: $0 <experiment_directory>"
    exit 1
fi

if [ ! -d "$EXPERIMENT_DIR" ]; then
    echo "Error: Experiment directory does not exist."
    echo "Usage: $0 <experiment_directory>"
    exit 1
fi

# Count BAM files

declare -a NORM_METHODS=("RPKM" "CPM" "BPM" "RPGC")
BAM_COUNT=$(find "${EXPERIMENT_DIR}/alignment" -maxdepth 1 -type f -name "*.sorted.bam" | wc -l)
echo "Found ${BAM_COUNT} BAM files"
TOTAL_JOBS=$((BAM_COUNT * ${#NORM_METHODS[@]}))
echo "Found ${TOTAL_JOBS} jobs to run"

if [ $BAM_COUNT -eq 0 ]; then
    echo "Error: No BAM files found in ${EXPERIMENT_DIR}/alignment"
    exit 1
fi

# Submit job
sbatch --array=1-${TOTAL_JOBS}%16 run_bamcoverage_normalizations.sh "$EXPERIMENT_DIR"
