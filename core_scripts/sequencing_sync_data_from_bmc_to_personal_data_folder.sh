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

# Ensure script runs in slurm.
CLUSTER_NAME="luria"
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
FILETYPE_TO_SYNC="*.fastq"
RSYNC_OPTS=(-av --include='*/' --include="$FILETYPE_TO_SYNC" --exclude='*')

#============================== 
# Setup and preprocessing
#============================== 
EXPERIMENT_DIR="$HOME/data/${EXPERIMENT_ID}"
DESTINATION_DIR="$EXPERIMENT_DIR/fastq/"
SOURCE_DIR="/net/bmc-pub17/data/bmc/public/Bell/${EXPERIMENT_ID}/"

#============================== 
# Error handling
#============================== 

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
