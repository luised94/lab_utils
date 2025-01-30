#!/bin/bash
# Script: submit_macs2_peak_calling
# Purpose: Submits run_macs2_peak_calling via slurm
# Usage: sbatch --array=1-N%16 run_bamcoverage_array.sbatch <experiment_directory>
# Dependencies: run_macs2_peak_calling, macs2, miniforge, slurm
# Date: 2024-12-26
# Usage: $./submit_macs2_peak_calling.sh $HOME/data/<experiment_directory>

# Immediately exit if not on luria.
if [[ "$(hostname)" != "luria" ]]; then
    echo "Error: This script must be run on luria cluster"
    exit 1
fi

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


# Count bam files
BAM_COUNT=$(find "${EXPERIMENT_DIR}/alignment" -maxdepth 1 -type f -name "*.bam" | wc -l)
echo "Found ${BAM_COUNT} bam files"

if [ $BAM_COUNT -eq 0 ]; then
    echo "Error: No bam files found in ${EXPERIMENT_DIR}/alignment"
    exit 1
fi

# Format file listing with columns and headers
echo -e "\nBAM files found:"
echo "----------------"
#mapfile -t unique_files < <(find "${BAM_DIRECTORY}" -maxdepth 1 -type f -name "*_sorted.bam" -exec basename {} \;)
#printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
#printf '%s\n' "${unique_files[@]}" | column -c "${COLUMNS:-$(tput cols)}"
#printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
find "${EXPERIMENT_DIR}/alignment" -maxdepth 1 -type f -name "*.bam" -exec basename {} \; | \
    pr -3 -t -w 100 | \
    column -t
echo "----------------"
echo -e "\nWill submit array job with following parameters:"
echo "Array size: 1-${BAM_COUNT}"
echo "Max simultaneous jobs: 16"
echo "Script: run_macs2_peak_calling.sbatch"
echo "Working directory: ${EXPERIMENT_DIR}"

read -p "Proceed with job submission? (y/n): " confirm
if [[ ! $confirm =~ ^[Yy]$ ]]; then
    echo "Job submission cancelled"
    exit 0
fi

# Submit job
sbatch --array=1-${BAM_COUNT}%16 "$HOME/lab_utils/core_scripts/run_macs2_peak_calling.sbatch" "$EXPERIMENT_DIR"
