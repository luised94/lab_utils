#!/bin/bash
################################################################################
# Submit chip processing script.
# Author: Luis | Date: 2025-10-20 | Version: 2.0.0
################################################################################
# PURPOSE:
#   Determine the read type (single vs paired), show job information and confirm with user to submit the CHIP processing pipeline.
# USAGE:
#   From the command line
#   $ ngs_submit_chip_processing.sh <EXPERIMENT_ID>
# DEPENDENCIES:
#   bash 4.2, bowtie2, fastp, bedtools, samtools
#   Assumes fastq files are in the fastq directory and everything is removed. (cleanup script.)
#   Assumes manifest file has been generated.
#   Processing script: ngs_run_chip_processing.sh
# OUTPUTS:
#   Sorted and indexed bam files, bigwig files for visualization.
################################################################################
# ============================================================================
# USAGE AND REQUIREMENTS
# ============================================================================
show_usage() {
    cat << EOF
Usage: $0 <experiment_directory>

Submit SLURM array job to process consolidated FASTQ files via the CHIP processing pipeline.

Arguments:
    EXPERIMENT_ID    Full path to experiment directory
                           Example: \$HOME/data/250930Bel
Options:
  -h, --help        Show this help message

Requirements:
    - Must run on luria cluster
    - Must be in git repository (lab_utils)
    - Consolidated FASTQ files must exist
    - Manifest file must exist: documentation/consolidated_reads_manifest.tsv

The script will:
  1. Read manifest to determine samples and read type
  2. Validate all FASTQ files exist
  3. Show processing summary
  4. Confirm with user.
  5. Submit SLURM array job

EOF
  exit 0
}

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

#==============================
# Argument error handling
#==============================
# Check for arguments
echo "Handling arguments..."
MIN_NUMBER_OF_ARGS=1
MAX_NUMBER_OF_ARGS=1
EXPECTED_EXPERIMENT_ID_PATTERN=^[0-9]{6}Bel$ # Do not quote regular expression.
# Check for help flag
if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
  show_usage
fi

# Verify number of arguments
if [[ $# -lt $MIN_NUMBER_OF_ARGS ]]; then
    echo "Error: Missing required argument EXPERIMENT_ID." >&2
    show_usage
    exit 1
fi

if [[ $# -gt $MAX_NUMBER_OF_ARGS ]]; then
    echo "Error: Too many arguments provided." >&2
    show_usage
    exit 1
fi

# Handle first argument: Remove trailing slash and validate pattern
EXPERIMENT_ID=${1%/} # Ensure argument does not have trailing slashes.
echo "Running error handling..."
if [[ ! $EXPERIMENT_ID =~ $EXPECTED_EXPERIMENT_ID_PATTERN ]]; then
  echo "Error: EXPERIMENT_ID does not match expected pattern." >&2
  echo "Please adjust EXPERIMENT_ID accordingly." >&2
  echo "EXPERIMENT ID PATTERN: $EXPECTED_EXPERIMENT_ID_PATTERN" >&2
  echo "EXPERIMENT_ID: $EXPERIMENT_ID" >&2
  exit 1
fi

# Handle second argument: Set dry-run mode.
#DRY_RUN=true
#if [[ $# -eq $MAX_NUMBER_OF_ARGS ]]; then
#    if [[ "$2" != "--active-run" ]]; then
#        echo "Error: Unknown option '$2'" >&2
#        echo "Use --active-run to perform actual sync." >&2
#        exit 1
#    fi
#    DRY_RUN=false
#fi

#==============================
# Configuration
#==============================
echo "Setting configuration..."
# Manifest output configuration
MANIFEST_FILENAME="consolidated_reads_manifest.tsv"
# SLURM job configuration
MAX_SIMULTANEOUS_JOBS=16
REPO_DIRECTORY_PATH=$(git rev-parse --show-toplevel)
PROCESSING_SCRIPT_TO_SUBMIT="${REPO_DIRECTORY_PATH}/core_scripts/ngs_run_chip_processing.sbatch"
NGS_FILE_PATTERNS=("consolidated*.fastq" "*.bam" "*.bw")

#==============================
# Setup and preprocessing
#==============================
EXPERIMENT_DIR="$HOME/data/${EXPERIMENT_ID}"
FASTQ_DIRECTORY="$EXPERIMENT_DIR/fastq/"
DOCUMENTATION_DIR="$(dirname "$FASTQ_DIRECTORY")/documentation"
JOB_LOG="${DOCUMENTATION_DIR}/experiment_job_info.md"
MANIFEST_FILEPATH="$DOCUMENTATION_DIR/$MANIFEST_FILENAME"
#FASTQ_FILE_PATTERN="consolidated*.fastq"

mkdir -p "$DOCUMENTATION_DIR" "$FASTQ_DIRECTORY"
touch "$JOB_LOG"

#==============================
# Error handling
#==============================
echo -e "\n=== CHIP Processing Submission Parameters ==="
echo "Working directory: ${EXPERIMENT_DIR}"
echo "Max simultaneous jobs: $MAX_SIMULTANEOUS_JOBS"
echo "Script: $PROCESSING_SCRIPT_TO_SUBMIT"

if [[ ! -d "$FASTQ_DIRECTORY" ]]; then
  echo "Error: FASTQ_DIRECTORY does not exist. Please verify experiment id." >&2
  echo "FASTQ_DIRECTORY: $FASTQ_DIRECTORY" >&2
  echo "Run $0 -h for additional help." >&2
  exit 1

fi

if [[ ! -d "$DOCUMENTATION_DIR" ]]; then
  echo "Error: DOCUMENTATION_DIR does not exist. Please verify experiment id." >&2
  echo "DOCUMENTATION_DIR: $DOCUMENTATION_DIR" >&2
  echo "Run $0 -h for additional help." >&2
  exit 1

fi

if [[ ! -f "$PROCESSING_SCRIPT_TO_SUBMIT" ]]; then
    echo "Error: PROCESSING_SCRIPT_TO_SUBMIT does not exist."
    echo "Value: $PROCESSING_SCRIPT_TO_SUBMIT"
    exit 1

fi

if [[ ! -f "$MANIFEST_FILEPATH" ]]; then
    echo "Error: Manifest file not found: $MANIFEST_FILEPATH"
    echo "The manifest should be created by the consolidation script."
    echo "Run consolidation first: consolidate_fastq_by_id.sh"
    exit 1
fi
# ============================================================================
# READ AND VALIDATE MANIFEST
# ============================================================================
# Count samples (skip header row)
TOTAL_SAMPLES=$(tail -n +2 "$MANIFEST_FILEPATH" | wc -l)

if [[ "$TOTAL_SAMPLES" -eq 0 ]]; then
    echo "Error: No samples found in manifest"
    echo "Manifest file: $MANIFEST_FILEPATH"
    exit 1
fi

# Detect read type from manifest structure
# Read first data row to check if read2_path column exists and has value
first_data_row=$(sed -n '2p' "$MANIFEST_FILEPATH")
IFS=$'\t' read -r sample_id read1_path read2_path <<< "$first_data_row"

if [[ -z "$read2_path" ]]; then
    read_type="single-end"
    echo "Detected read type: single-end"

else
    read_type="paired-end"
    echo "Detected read type: paired-end"

fi

#==============================
# Output experiment information for user confirmation
#==============================
# Count different file patterns.
FASTQ_COUNT=0

for pattern in "${NGS_FILE_PATTERNS[@]}"; do
    count=$(find "$EXPERIMENT_DIR" -type f -name "$pattern" 2>/dev/null | wc -l)

    if [[ "$pattern" == *"fastq"* ]]; then
        FASTQ_COUNT=$count
        echo "Found $count FASTQ files matching pattern: $pattern"
    else
        # Derive type from extension for message
        ext="${pattern#*.}"
        echo "Found $count $ext files matching pattern: $pattern"
    fi
done

# Script processes starting from fastq files.
if [[ "$FASTQ_COUNT" -eq 0 ]]; then
    echo "Error: No fastq files found in ${FASTQ_DIR}"
    exit 1

fi


echo "Array size: 1-${FASTQ_COUNT}"
read -rp "Proceed with job submission? (y/n): " confirm
if [[ ! $confirm =~ ^[Yy]$ ]]; then
    echo "Job submission cancelled"
    exit 0

fi

# --- Capture job description ---
#description=""
#while [[ -z "$description" ]]; do
#    read -rp "Enter job purpose (e.g., 'Filter reads for sample X'): " description
#    if [[ -z "$description" ]]; then
#        echo "Error: Description cannot be empty" >&2
#    fi
#
#done
#
## --- Submit job and capture output ---
#job_submit_output=$(sbatch --array=1-"${FASTQ_COUNT}%$MAX_SIMULTANEOUS_JOBS" \
#    "$PROCESSING_SCRIPT_TO_SUBMIT" \
#    "$EXPERIMENT_DIR")
#
## Extract job ID (works with both "Submitted batch job 12345" and "12345")
#job_id=$(echo "$job_submit_output" | grep -oE '[0-9]+$')
#
## --- Log job details ---
#{
#    echo "# $job_id"
#    echo "- Submission time: $(date --iso-8601=seconds)"
#    echo "- Cluster: $(hostname)"
#    # Git metadata
#    (
#        # Ensure the git commands are executed inside the repository.
#        cd "$HOME/lab_utils" || exit
#        echo "- Git commit: $(git rev-parse --short HEAD 2>/dev/null || echo 'unknown')"
#        echo "- Git branch: $(git symbolic-ref --short HEAD 2>/dev/null || echo 'detached')"
#        echo "- Git status: $(git status --porcelain 2>/dev/null | wc -l) uncommitted changes"
#    )
#    echo "- Experiment dir: $EXPERIMENT_DIR"
#    echo "- Command ran: $0"
#    echo "- sbatch command: sbatch --array=1-${FASTQ_COUNT}%$MAX_SIMULTANEOUS_JOBS $PROCESSING_SCRIPT_TO_SUBMIT $EXPERIMENT_DIR"
#    echo "- FASTQ files processed: $FASTQ_COUNT"
#    echo "- Description: $description"
#    echo "- Logs: {{fill out comments}}"
#    echo ""
#} >> "$JOB_LOG"
#
#echo "Job $job_id submitted successfully. Details logged to $JOB_LOG"

# Count fastq files
#mapfile -t PREFILTERED_FASTQ_FILENAMES < <(find "${FASTQ_DIR}" -maxdepth 1 -type f -name "$FASTQ_FILE_PATTERN" -exec basename {} \;)
#FASTQ_COUNT=${#PREFILTERED_FASTQ_FILENAMES[@]}
#echo "Found ${FASTQ_COUNT} fastq files"
# Format file listing with columns and headers
#echo -e "\nFASTQ files found:"
#printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
#printf '%s\n' "${PREFILTERED_FASTQ_FILENAMES[@]}" | column -c "${COLUMNS:-$(tput cols)}"
#printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
