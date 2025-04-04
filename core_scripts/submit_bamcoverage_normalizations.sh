#!/bin/bash
# submit_bamcoverage_normalizations.sh
# Usage: $./submit_bamcoverage_normalizations.sh $HOME/data/<experiment_directory>
#
# Immediately exit if not on luria.
if [[ "$(hostname)" != "luria" ]]; then
    echo "Error: This script must be run on luria cluster"
    exit 1
fi

EXPERIMENT_DIR="$1"
if [ -z "$EXPERIMENT_DIR" ]; then
    echo "Error: Experiment directory not provided"
    echo "Usage: $0 <experiment_directory>"
    echo "Ensure experiment directory is full path: $HOME/data/<experiment_directory"
    exit 1
fi

if [ ! -d "$EXPERIMENT_DIR" ]; then
    echo "Error: Experiment directory does not exist."
    echo "Usage: $0 <experiment_directory>"
    exit 1
fi

BAM_DIRECTORY="${EXPERIMENT_DIR}/alignment"

# Count BAM files
BAM_COUNT=$(find "${BAM_DIRECTORY}" -maxdepth 1 -type f -name "*_sorted.bam" | wc -l)
echo "Found ${BAM_COUNT} BAM files"

declare -a NORM_METHODS=("RPKM" "CPM" "BPM" "RPGC")
TOTAL_JOBS=$((BAM_COUNT * ${#NORM_METHODS[@]}))
echo "Found ${TOTAL_JOBS} jobs to run"
if [ $BAM_COUNT -eq 0 ]; then
    echo "Error: No BAM files found in ${BAM_DIRECTORY}"
    exit 1
fi

mapfile -t unique_files < <(find "${BAM_DIRECTORY}" -maxdepth 1 -type f -name "processed*_sorted.bam" -exec basename {} \;)

# Format file listing with columns and headers
echo -e "\nBAM files found:"
printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
printf '%s\n' "${unique_files[@]}" | column -c "${COLUMNS:-$(tput cols)}"
printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -

echo -e "\nWill submit array job with following parameters:"
echo "Array size: 1-${BAM_COUNT}"
echo "Max simultaneous jobs: 16"
echo "Script: run_bamcoverage_normalizations.sbatch"
echo "Working directory: ${EXPERIMENT_DIR}"

read -p "Proceed with job submission? (y/n): " confirm
if [[ ! $confirm =~ ^[Yy]$ ]]; then
    echo "Job submission cancelled"
    exit 0
fi

# Submit job
sbatch --array=1-${TOTAL_JOBS}%16 "$HOME/lab_utils/core_scripts/run_bamcoverage_normalizations.sbatch" "$EXPERIMENT_DIR"
