#!/bin/bash
# bmc_submit_chip_processing.sh
# Dependencies:
# - Single end: Cleaned up and consolidated fastq files.
# - Paired end: Validate properly paired reads
# Usage: $./bmc_submit_chip_processing.sh ./data/<experiment_directory>

# Immediately exit if not on luria.
if [[ "$(hostname)" != "luria" ]]; then
    echo "Error: This script must be run on luria cluster"
    exit 1

fi

if ! git rev-parse --git-dir > /dev/null 2>&1; then
    echo "Error: Not in a git repository"
    echo "Move to lab_utils directory."
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

FASTQ_FILE_PATTERN="consolidated*.fastq"
FASTQ_DIR="${EXPERIMENT_DIR}/fastq"
DOCUMENTATION_DIR="${EXPERIMENT_DIR}/documentation"
JOB_LOG="${DOCUMENTATION_DIR}/experiment_job_info.md"
REPO_DIRECTORY_PATH=$(git rev-parse --show-toplevel)
SCRIPT_TO_SUBMIT="${REPO_DIRECTORY_PATH}/lab_utils/core_scripts/bmc_run_chip_processing.sbatch"
touch "$JOB_LOG"
mkdir -p "$DOCUMENTATION_DIR" "$FASTQ_DIR"

# Count fastq files
mapfile -t PREFILTERED_FASTQ_FILENAMES < <(find "${FASTQ_DIR}" -maxdepth 1 -type f -name "$FASTQ_FILE_PATTERN" -exec basename {} \;)
FASTQ_COUNT=${#PREFILTERED_FASTQ_FILENAMES[@]}
echo "Found ${FASTQ_COUNT} fastq files"

if [[ "$FASTQ_COUNT" -eq 0 ]]; then
    echo "Error: No fastq files found in ${FASTQ_DIR}"
    exit 1

fi

if [[ ! -f "$SCRIPT_TO_SUBMIT" ]]; then
    echo "Error: SCRIPT_TO_SUBMIT does not exist."
    echo "Value: $SCRIPT_TO_SUBMIT"
    exit 1

fi

# Format file listing with columns and headers
echo -e "\nFASTQ files found:"
printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
printf '%s\n' "${PREFILTERED_FASTQ_FILENAMES[@]}" | column -c "${COLUMNS:-$(tput cols)}"
printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
echo -e "\nWill submit array job with following parameters:"
echo "Array size: 1-${FASTQ_COUNT}"
echo "Max simultaneous jobs: 16"
echo "Script: $SCRIPT_TO_SUBMIT"
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
    "$SCRIPT_TO_SUBMIT" \
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
    echo "- sbatch command: sbatch --array=1-${FASTQ_COUNT}%16 $SCRIPT_TO_SUBMIT $EXPERIMENT_DIR"
    echo "- FASTQ files processed: $FASTQ_COUNT"
    echo "- Description: $description"
    echo "- Logs: {{fill out comments}}"
    echo ""
} >> "$JOB_LOG"

echo "Job $job_id submitted successfully. Details logged to $JOB_LOG"
