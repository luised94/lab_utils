#!/bin/bash
# submit_alignment.sh
# Dependencies: Filtered files by fastp
# Usage: $./submit_alignment.sh $HOME/data/<experiment_directory>
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
    exit 1
fi

if [ ! -d "$EXPERIMENT_DIR" ]; then
    echo "Error: Experiment directory does not exist."
    echo "Usage: $0 <experiment_directory>"
    exit 1
fi

FASTQ_DIR="${EXPERIMENT_DIR}/fastq"

# Count fastq files
FASTQ_COUNT=$(find "${FASTQ_DIR}" -maxdepth 1 -type f -name "processed*.fastq" | wc -l)
echo "Found ${FASTQ_COUNT} fastq files"

if [ $FASTQ_COUNT -eq 0 ]; then
    echo "Error: No fastq files found in $FASTQ_DIR"
    exit 1
fi

mapfile -t unique_files < <(find "${FASTQ_DIR}" -maxdepth 1 -type f -name "processed*.fastq" -exec basename {} \;)

# Format file listing with columns and headers
echo -e "\nFASTQ files found:"
printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
printf '%s\n' "${unique_files[@]}" | column -c "${COLUMNS:-$(tput cols)}"
printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
#find "${EXPERIMENT_DIR}/fastq" -maxdepth 1 -type f -name "processed*.fastq" -exec basename {} \; | \
#    pr -3 -t -w 100 | \
#    column -t

echo -e "\nWill submit array job with following parameters:"
echo "Array size: 1-${FASTQ_COUNT}"
echo "Max simultaneous jobs: 16"
echo "Script: run_bowtie2_array_alignment.sbatch"
echo "Working directory: ${EXPERIMENT_DIR}"

read -p "Proceed with job submission? (y/n): " confirm
if [[ ! $confirm =~ ^[Yy]$ ]]; then
    echo "Job submission cancelled"
    exit 0
fi

# Submit job
sbatch --array=1-${FASTQ_COUNT}%16 "$HOME/lab_utils/core_scripts/run_bowtie2_array_alignment.sbatch" "$EXPERIMENT_DIR"
