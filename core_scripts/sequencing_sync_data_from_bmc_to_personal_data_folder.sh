#!/bin/bash
################################################################################
# Sync the data from the BMC's long term storage to personal folder for use
# Author: Luis | Date: 2025-10-20 | Version: 2.0.0
################################################################################
# PURPOSE: 
#   For an experiment id, sync the fastq files from the long term storage to the personal data folder for analysis
# USAGE:
#   From the command line and while working on the cluster.
#   $ ./sequencing_sync_data_from_bmc_to_personal_data_folder.sh <EXPERIMENT_ID>
# DEPENDENCIES:
#   bash 4.2
#   slurm environment (srun)
# OUTPUTS:
#   EXPERIMENT DIR with fastq files with the directory structure 
################################################################################
#============================== 
# Usage and help
#============================== 
show_usage() {
  cat << EOF
Usage: $(basename "$0") <fastq_directory> [-v]

Description:
  Rsync the data (not the directory)from the BMC long term storage to local data folder for analysis.
  Script runs in dry-run by default. Provide --active-run option to execute.

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

if [[ "$2" != "--active-run" ]]; then
  echo "Error: Invalid second argument: '$2'" >&2
  echo "Expected: --active-run (or omit for dry-run)" >&2
  echo "Run '$0 -h' for usage." >&2
  exit 1
else
  DRY_RUN=false
fi

# Ensure script runs in slurm.
if [[ ! "$(hostname)" == "$CLUSTER_NAME" ]]; then
    echo "Error: This script should be run on the cluster." >&2
    echo "Current cluster setting: $CLUSTER_NAME" >&2
    echo "Please ssh into cluster." >&2
    exit 1

fi

#============================== 
# Configuration
#============================== 
echo "Setting configuration..."
DRY_RUN=true
CLUSTER_NAME="luria"
FILETYPE_TO_SYNC="*.fastq"
EXPECTED_EXPERIMENT_ID_PATTERN=^[0-9]{8}Bel$ # Do not quote regular expression.
RSYNC_OPTS=(-av --include='*/' --include="$FILETYPE_TO_SYNC" --exclude='*')

#============================== 
# Setup and preprocessing
#============================== 
# Ensure argument does not have trailing slashes.
EXPERIMENT_ID=${1%/}
EXPERIMENT_DIR="$HOME/data/${EXPERIMENT_ID}"
DESTINATION_DIR="$EXPERIMENT_DIR/fastq/"
SOURCE_DIR="/net/bmc-pub17/data/bmc/public/Bell/${EXPERIMENT_ID}/"

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

if [[ ! -d "$SOURCE_DIR" ]]; then
  echo "Error: SOURCE_DIR does not exist. Please adjust parameter." >&2
  echo "Verify that proper bmc-pub parameter or experiment id." >&2
  echo "SOURCE_DIR: $SOURCE_DIR" >&2
  echo "Please adjust EXPERIMENT_ID accordingly." >&2
  exit 1

fi

mapfile -t files_to_sync < <(find "$SOURCE_DIR" -type f -name "$FILETYPE_TO_SYNC" )
if [[ ${#files_to_sync[@]} -eq 0 ]]; then
  echo "No files to sync found." >&2
  echo "FILETYPE_TO_SYNC: $FILETYPE_TO_SYNC" >&2
  echo "SOURCE_DIR: $SOURCE_DIR" >&2
  exit 1

fi

mkdir -p "$DESTINATION_DIR"

#============================== 
# Main logic
#============================== 
echo "Executing rsync command..."
if [[ "$DRY_RUN" == true ]]; then
    RSYNC_OPTS=(-n "${RSYNC_OPTS[@]}")  # add dry-run flag
    echo ">>> DRY RUN MODE (no files will be copied)" >&2
fi

srun rsync "${RSYNC_OPTS[@]}" "$SOURCE_DIR" "$DESTINATION_DIR"
