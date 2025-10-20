#!/bin/bash
################################################################################
# Verify paired end fastq files and output manifest
# Author: Luis | Date: 2025-10-20 | Version: 1.0.0
################################################################################
# PURPOSE:
#   Remove unused files from rsync download and moves all fastq files out of subdirectories.
# USAGE:
#   From the command line
#   $ srun ./sequencing_verify_fastq_pairs.sh /path/to/fastq/directory
# DEPENDENCIES:
#   bash 4.2
# OUTPUTS:
#   No outputs produced by file or script.
#   Scripts (such as setup_bmc_experiment.R) source this script.
################################################################################
# Dependencies: Requires fastq files downloaded to the appropriate experiment fastq directory. Follow instructions in email with subject * Data Ready from the BMC personnel.
# srun rsync -av /net/bmc-pub17/data/bmc/public/Bell/${PROJECT_ID}/ ~/data/${PROJECT_ID}/fastq/
# Usage: Run relevant code lines manually. Not worth it to automate currently.

# Strict error handling
set -euo pipefail
trap 'echo "Error on line $LINENO"' ERR

# Log file in /tmp for operations tracking
log_file="/tmp/cleanup_$(date +%Y%m%d_%H%M%S).log"
exec 1> >(tee -a "$log_file")
exec 2>&1

echo "Starting cleanup operation at $(date)"

# Store current directory
current_dir=$(pwd)
echo "Working directory: $current_dir"

# First count existing fastq files for verification
initial_fastq_count=$(find . -type f -name "*.fastq" | wc -l)
echo "Found $initial_fastq_count FASTQ files initially"
#[[  $initial_fastq_count -eq 0 ]] && { echo "No fastq files in directory."; return; }
[[  $initial_fastq_count -eq 0 ]] && { echo "No fastq files in directory."; }

# You can perform a dry-run first by removing the -delete option from the command.
# Use | wc -l to count the files as well.
# Remove unmapped files first
echo "Removing unmapped files..."
find . -type f -name "*unmapped*" -delete

# Remove all non-fastq files
echo "Removing non-FASTQ files..."
find . -type f ! -name "*.fastq" -delete

# Move all fastq files to current directory
echo "Moving FASTQ files to current directory..."
find . -type f -name "*.fastq" -exec mv {} . \;

# Remove empty directories
echo "Removing empty directories..."
find . -type d -empty -delete

# Verify final state
# Should be double if you used Aviti. Two for each lane.
final_fastq_count=$(find . -maxdepth 1 -type f -name "*.fastq" | wc -l)
echo "Final FASTQ count in current directory: $final_fastq_count"

# This comparison does not work. Initially there are two files per sample and the additional fastq for each lane.
#if [ "$initial_fastq_count" -ne "$final_fastq_count" ]; then
#    echo "ERROR: FASTQ file count mismatch! Initial: $initial_fastq_count, Final: $final_fastq_count"
#    exit 1
#fi

echo "Operation completed successfully at $(date)"
echo "Log file: $log_file"
