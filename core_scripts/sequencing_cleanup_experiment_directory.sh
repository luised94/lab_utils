#!/bin/bash
################################################################################
# Verify paired end fastq files and output manifest
# Author: Luis | Date: 2025-10-20 | Version: 1.0.0
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
if [[ $# -eq 0 ]]; then
  echo "Error: No argument provided." >&2
  show_usage
fi

# Check for help flag
if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
  show_usage
fi

if [[ $# -gt 2 ]]; then
    echo "Error: Too many arguments provided." >&2
    echo "Run '$0 -h' for usage." >&2
    exit 1
fi

DRY_RUN=true
if [[ "$2" != "--active-run" ]]; then
  echo "Error: Invalid second argument: '$2'" >&2
  echo "Expected: --active-run (or omit for dry-run)" >&2
  echo "Run '$0 -h' for usage." >&2
  exit 1
else
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
EXPECTED_EXPERIMENT_ID_PATTERN=^[0-9]{8}Bel$ # Do not quote regular expression.

#============================== 
# Setup and preprocessing
#============================== 
# Ensure argument does not have trailing slashes.
EXPERIMENT_ID=${1%/}
EXPERIMENT_DIR="$HOME/data/${EXPERIMENT_ID}"
FASTQ_DIR="$EXPERIMENT_DIR/fastq/"

#============================== 
# Error handling
#============================== 
echo "Running error handling..."
if [[ ! $EXPERIMENT_ID =~ $EXPECTED_EXPERIMENT_ID_PATTERN ]]; then
  echo "Error: EXPERIMENT_ID does not match expected pattern." >&2
  echo "Please adjust EXPERIMENT_ID accordingly." >&2
  echo "EXPERIMENT ID PATTERN: $EXPECTED_EXPERIMENT_ID_PATTERN" >&2
  echo "EXPERIMENT_ID: $EXPERIMENT_ID" >&2
  exit 1
fi

if [[ ! -d "$FASTQ_DIR" ]]; then
  echo "Error: FASTQ_DIR does not exist. Please verify experiment id." >&2
  echo "FASTQ_DIR: $FASTQ_DIR" >&2
  echo "Run $0 -h for additional help." >&2
  exit 1

fi

# ############################################
# Main logic
# ############################################

mapfile -t files_to_move < <(
  find "$FASTQ_DIR" \
       -type f \
       -name "$FILETYPE_TO_KEEP" \
       -not -name "$EXCLUDING_PATTERN"
  )

mapfile -t files_to_remove < <(
  find "$FASTQ_DIR" \
       ! -type f \
       -name "$FILETYPE_TO_KEEP" \
       -not -name "$EXCLUDING_PATTERN"

  )

FASTQ_FILES_COUNT=${#files_to_move[@]}
if [[ $FASTQ_FILES_COUNT -eq 0 ]]; then
  echo "No fastq files found." >&2
  echo "FASTQ_DIR: $FASTQ_DIR" >&2
  exit 1

fi

#if [[ ${#files_to_remove[@]} -eq 0 ]]; then
#  echo "No files to remove found." >&2
#  echo "EXCLUDING_PATTERN: $EXCLUDING_PATTERN" >&2
#  exit 1
#
#fi

if [[ $DRY_RUN==true ]]; then
  echo "Executing dry-run..."
fi
# You can perform a dry-run first by removing the -delete option from the command.
# Use | wc -l to count the files as well.
# Remove unmapped files first
#echo "Removing unmapped files..."
#find . -type f -name "*unmapped*" -delete
#
## Remove all non-fastq files
#echo "Removing non-FASTQ files..."
#find . -type f ! -name "*.fastq" -delete
#
## Move all fastq files to current directory
#echo "Moving FASTQ files to current directory..."
#find . -type f -name "*.fastq" -exec mv {} . \;
#
## Remove empty directories
#echo "Removing empty directories..."
#find . -type d -empty -delete
