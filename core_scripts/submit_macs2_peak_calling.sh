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
    echo "Ensure experiment directory is full path: $HOME/data/<experiment_directory"
    exit 1
fi

if [ ! -d "$EXPERIMENT_DIR" ]; then
    echo "Error: Experiment directory does not exist."
    echo "Usage: $0 <experiment_directory>"
    exit 1
fi

BAM_DIRECTORY="${EXPERIMENT_DIR}/alignment"
DOCUMENTATION_DIR="${EXPERIMENT_DIR}/documentation"
JOB_LOG="${DOCUMENTATION_DIR}/slurm_job_info.md"
mkdir -p "$DOCUMENTATION_DIR" "$BAM_DIRECTORY"

# Count bam files
BAM_COUNT=$(find "${BAM_DIRECTORY}" -maxdepth 1 -type f -name "processed*_sorted.bam" | wc -l)
echo "Found ${BAM_COUNT} bam files"

if [ "$BAM_COUNT" -eq 0 ]; then
    echo "Error: No bam files found in $BAM_DIRECTORY"
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
echo "Script: run_macs2_peak_calling.sbatch"
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
job_submit_output=$(sbatch --array=1-"${BAM_COUNT}%16" \
    "$HOME/lab_utils/core_scripts/run_macs2_peak_calling.sbatch" \
    "$EXPERIMENT_DIR")

# Extract job ID (works with both "Submitted batch job 12345" and "12345")
job_id=$(echo "$job_submit_output" | grep -oE '[0-9]+$')

# --- Log job details ---
{
    echo "# $job_id"
    echo "- Submission time: $(date --iso-8601=seconds)"
    echo "- Cluster: $(hostname)"
    # Git metadata
    (
        # Ensure the git commands are executed inside the repository.
        cd "$HOME/lab_utils" || exit
        echo "- Git commit: $(git rev-parse --short HEAD 2>/dev/null || echo 'unknown')"
        echo "- Git branch: $(git symbolic-ref --short HEAD 2>/dev/null || echo 'detached')"
        echo "- Git status: $(git status --porcelain 2>/dev/null | wc -l) uncommitted changes"
    )
    echo "- Experiment dir: $EXPERIMENT_DIR"
    echo "- Command ran: $0"
    echo "- sbatch command: sbatch --array=1-${TOTAL_JOBS}%16 $HOME/lab_utils/core_scripts/run_macs2_peak_calling.sbatch $EXPERIMENT_DIR"
    echo "- Files processed: $BAM_COUNT"
    echo "- Description: $description"
    echo "- Logs: {{fill out comments}}"
    echo ""
} >> "$JOB_LOG"
echo "Job $job_id submitted successfully. Details logged to $JOB_LOG"
