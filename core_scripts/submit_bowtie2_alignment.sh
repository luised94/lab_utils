#!/bin/bash
# submit_alignment.sh
# Dependencies: Filtered files by fastp
# Usage: $./submit_bowtie2_alignment.sh $HOME/data/<experiment_directory>
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

FASTQ_DIR="${EXPERIMENT_DIR}/fastq"
DOCUMENTATION_DIR="${EXPERIMENT_DIR}/documentation"
JOB_LOG="${DOCUMENTATION_DIR}/slurm_job_info.md"
mkdir -p "$DOCUMENTATION_DIR" "$FASTQ_DIR"

# Count fastq files
FASTQ_COUNT=$(find "${FASTQ_DIR}" -maxdepth 1 -type f -name "processed*.fastq" | wc -l)
echo "Found ${FASTQ_COUNT} fastq files"

if [ "$FASTQ_COUNT" -eq 0 ]; then
    echo "Error: No fastq files found in $FASTQ_DIR"
    exit 1
fi

mapfile -t unique_files < <(find "${FASTQ_DIR}" -maxdepth 1 -type f -name "processed*.fastq" -exec basename {} \;)

# Format file listing with columns and headers
echo -e "\nFASTQ files found:"
printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
printf '%s\n' "${unique_files[@]}" | column -c "${COLUMNS:-$(tput cols)}"
printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
echo -e "\nWill submit array job with following parameters:"
echo "Array size: 1-${FASTQ_COUNT}"
echo "Max simultaneous jobs: 16"
echo "Script: run_bowtie2_array_alignment.sbatch"
echo "Working directory: ${EXPERIMENT_DIR}"

read -rp "Proceed with job submission? (y/n): " confirm
if [[ ! $confirm =~ ^[Yy]$ ]]; then
    echo "Job submission cancelled"
    exit 0
fi

# --- Capture job description ---
description=""
while [[ -z "$description" ]]; do
    read -rp "Enter job purpose (e.g., 'Filter reads for sample X'): " description
    if [[ -z "$description" ]]; then
        echo "Error: Description cannot be empty" >&2
    fi
done

# --- Submit job and capture output ---
job_submit_output=$(sbatch --array=1-"${FASTQ_COUNT}%16" \
    "$HOME/lab_utils/core_scripts/run_bowtie2_array_alignment.sbatch" \
    "$EXPERIMENT_DIR")

# Extract job ID (works with both "Submitted batch job 12345" and "12345")
job_id=$(echo "$job_submit_output" | grep -oE '[0-9]+$')

# --- Log job details ---
{
    echo "# $job_id"
    echo "- Submission time: $(date --iso-8601=seconds)"
    echo "- Cluster: $(hostname)"
    echo "- Experiment dir: $EXPERIMENT_DIR"
    echo "- Command ran: $0"
    echo "- sbatch command: sbatch --array=1-${FASTQ_COUNT}%16 $HOME/lab_utils/core_scripts/run_bowtie2_array_alignment.sbatch $EXPERIMENT_DIR"
    echo "- FASTQ files processed: $FASTQ_COUNT"
    echo "- Description: $description"
    echo ""
} >> "$JOB_LOG"
echo "Job $job_id submitted successfully. Details logged to $JOB_LOG"
