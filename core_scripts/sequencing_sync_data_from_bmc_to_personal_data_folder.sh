#!/bin/bash
################################################################################
# Sync the data from the BMC's long term storage to personal folder for use
# Author: Luis | Date: 2025-10-20 | Version: 1.0.0
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

CLUSTER_NAME="luria"
if [[ ! "$(hostname)" == "$CLUSTER_NAME" ]]; then
    echo "Error: This script should be run on the cluster." >&2
    echo "Current cluster setting: $CLUSTER_NAME" >&2
    echo "Please ssh into cluster." >&2
    exit 1

fi

# Configuration
FILETYPE_TO_SYNC="*.fastq"
EXPECTED_EXPERIMENT_ID_PATTERN=^[0-9]{8}Bel$ # Do not quote regular expression.

# Handle arguments
# Ensure argument does not have trailing slashes.
EXPERIMENT_ID=${1%/}
EXPERIMENT_DIR="$HOME/data/${EXPERIMENT_ID}"
DESTINATION_DIR="$EXPERIMENT_DIR/fastq/"
SOURCE_DIR="/net/bmc-pub17/data/bmc/public/Bell/${EXPERIMENT_ID}/"

if [[ ! $EXPERIMENT_ID =~ $EXPECTED_EXPERIMENT_ID_PATTERN ]]; then
  echo "Error: EXPERIMENT_ID does not match expected pattern." >&2
  echo "Please adjust EXPERIMENT_ID accordingly." >&2
  echo "EXPERIMENT ID PATTERN: $EXPECTED_EXPERIMENT_ID_PATTERN" >&2
  echo "EXPERIMENT_ID: $EXPERIMENT_ID" >&2
  exit 1
fi

mkdir -p "$DESTINATION_DIR"

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

# Make verbose
#echo "------- Current Settings -------"
#echo "No files to sync found." >&2
#echo "FILETYPE_TO_SYNC: $FILETYPE_TO_SYNC" >&2
#echo "SOURCE_DIR: $SOURCE_DIR" >&2
#
#echo "------- Current Settings -------"

srun rsync -nav \
           --include=*/ \
           --include="$FILETYPE_TO_SYNC" \
           --exclude="*" \
           "$SOURCE_DIR" \
           "$DESTINATION_DIR"
