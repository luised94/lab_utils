#!/bin/bash
################################################################################
# Verify paired end fastq files and output manifest
# Author: Luis | Date: 2025-10-20 | Version: 2.0.0
################################################################################
# PURPOSE:
#   Remove unused files from rsync download and moves all fastq files out of subdirectories.
# USAGE:
#   From the command line
#   $ srun sequencing_cleanup_experiment_directory.sh <EXPERIMENT_ID>
# DEPENDENCIES:
#   bash 4.2
#   Assumes fastqs were downloaded using rsync.
# OUTPUTS:
#   Fastq files in the ~/data/<EXPERIMENT_ID>/fastq
#   Removes other file types and unmapped files.
################################################################################
#============================== 
# Usage and help
#============================== 
show_usage() {
  cat << EOF
Usage: srun $(basename "$0") <fastq_directory> [-v]

Description:
  Remove unused non-fastq files from fastq directory
  Script runs in dry-run by default. Provide --active-run option to execute, after reviewing output messages.

Arguments:
  EXPERIMENT_ID    Experiment id of experiment (from BMC submission.)
                   (e.g., 250930Bel)

Options:
  -h, --help        Show this help message
  --active-run      Execute rsync command

Output:
  Directory with fastq files downloaded to the local directly in initial directory structure.

Example:
  $(basename "$0") 250930Bel
  $(basename "$0") 250930Bel --active-run
  $(basename "$0") -h

EOF
  exit 0
}

#============================== 
# Argument error handling
#============================== 
# Check for arguments
echo "Handling arguments..."
MIN_NUMBER_OF_ARGS=1
MAX_NUMBER_OF_ARGS=2
EXPECTED_EXPERIMENT_ID_PATTERN=^[0-9]{8}Bel$ # Do not quote regular expression.
# Check for help flag
if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
  show_usage
fi

# Verify number of arguments
if [[ $# -lt $MIN_NUMBER_OF_ARGS ]]; then
    echo "Error: Missing required argument EXPERIMENT_ID." >&2
    show_usage
    exit 1
elif [[ $# -gt $MAX_NUMBER_OF_ARGS ]]; then
    echo "Error: Too many arguments provided ($#)." >&2
    show_usage
    exit 1
fi

# Handle first argument: Remove trailing slash and validate pattern
EXPERIMENT_ID=${1%/}
echo "Running error handling..."
if [[ ! $EXPERIMENT_ID =~ $EXPECTED_EXPERIMENT_ID_PATTERN ]]; then
  echo "Error: EXPERIMENT_ID does not match expected pattern." >&2
  echo "Please adjust EXPERIMENT_ID accordingly." >&2
  echo "EXPERIMENT ID PATTERN: $EXPECTED_EXPERIMENT_ID_PATTERN" >&2
  echo "EXPERIMENT_ID: $EXPERIMENT_ID" >&2
  exit 1
fi

# Handle second argument: Set dry-run mode.
DRY_RUN=true
if [[ $# -eq $MAX_NUMBER_OF_ARGS ]]; then
    if [[ "$2" != "--active-run" ]]; then
        echo "Error: Unknown option '$2'" >&2
        echo "Use --active-run to perform actual sync." >&2
        exit 1
    fi
    DRY_RUN=false
fi

# Ensure script is run inside a Slurm allocation
if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    cat >&2 <<EOF
Error: This script must be run within a Slurm job.

To run interactively:
    srun $0 <EXPERIMENT_ID> [--active-run]

To submit as a batch job:
    echo "$0 <EXPERIMENT_ID>" [--active-run] | sbatch

EOF
    exit 1
fi

#============================== 
# Configuration
#============================== 
echo "Setting configuration..."
FILETYPE_TO_KEEP="*.fastq"
EXCLUDING_PATTERN="*unmapped*"

#============================== 
# Setup and preprocessing
#============================== 
EXPERIMENT_DIR="$HOME/data/${EXPERIMENT_ID}"
FASTQ_DIRECTORY="$EXPERIMENT_DIR/fastq/"
mkdir -p "$FASTQ_DIRECTORY"

#============================== 
# Error handling
#============================== 
if [[ ! -d "$FASTQ_DIRECTORY" ]]; then
  echo "Error: FASTQ_DIRECTORY does not exist. Please verify experiment id." >&2
  echo "FASTQ_DIRECTORY: $FASTQ_DIRECTORY" >&2
  echo "Run $0 -h for additional help." >&2
  exit 1

fi

# ############################################
# Main logic
# ############################################
FASTQ_FILES_COUNT=$(find "$FASTQ_DIRECTORY" -type f -name "$FILETYPE_TO_KEEP" ! -name "$EXCLUDING_PATTERN" 2>/dev/null | wc -l)
if [[ $FASTQ_FILES_COUNT -eq 0 ]]; then
  echo "No fastq files found." >&2
  echo "FASTQ_DIRECTORY: $FASTQ_DIRECTORY" >&2
  exit 1

fi

if [[ "$DRY_RUN" == true ]]; then
  echo ">>> DRY RUN MODE <<<"
  echo "Files that would be MOVED to $FASTQ_DIRECTORY:"
  # Files to move: valid FASTQ files that are NOT in the top level
  mapfile -t to_move < <(
      find "$FASTQ_DIRECTORY" -mindepth 2 -type f -name "$FILETYPE_TO_KEEP" ! -name "$EXCLUDING_PATTERN" 2>/dev/null
  )

  # Files to move (valid FASTQ)
  if [[ ${#to_move[@]} -eq 0 ]]; then
    echo "No FASTQ files to move (all already in top level or none found)."

  else
    echo "Files that would be moved to: $FASTQ_DIRECTORY"
    printf "%s\n" "${to_move[@]}"

  fi

  # Files to delete (non-FASTQ OR unmapped FASTQ)
  mapfile -t to_delete < <(find "$FASTQ_DIRECTORY" -type f \( ! -name "$FILETYPE_TO_KEEP" -o -name "$EXCLUDING_PATTERN" \) 2>/dev/null)
  if [[ ${#to_delete[@]} -eq 0 ]]; then
    echo "No unwanted files to delete."

  else
    echo "Files that would be deleted:"
    printf "%s\n" "${to_delete[@]}"

  fi

  # Empty directories
  mapfile -t empty_dirs < <(find "$FASTQ_DIRECTORY" -type d -empty 2>/dev/null)
  if [[ ${#empty_dirs[@]} -eq 0 ]]; then
    echo "No empty directories to remove."

  else
    echo "Empty directories that would be removed:"
    printf "%s\n" "${empty_dirs[@]}"

  fi

else
  echo "Cleaning up directory..."
  # 1. Move only files from subdirectories
  find "$FASTQ_DIRECTORY" -mindepth 2 -type f -name "$FILETYPE_TO_KEEP" ! -name "$EXCLUDING_PATTERN" \
     -exec mv {} "$FASTQ_DIRECTORY" \;

  # 2. Delete everything else (non-fastq OR unmapped fastq)
  find "$FASTQ_DIRECTORY" -type f \( ! -name "$FILETYPE_TO_KEEP" -o -name "$EXCLUDING_PATTERN" \) \
       -delete

  # 3. Clean empty dirs
  find "$FASTQ_DIRECTORY" -type d -empty -delete

  echo "Cleanup complete."

fi
echo "Script done..."
