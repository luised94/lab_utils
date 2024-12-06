#!/bin/bash

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
final_fastq_count=$(find . -maxdepth 1 -type f -name "*.fastq" | wc -l)
echo "Final FASTQ count in current directory: $final_fastq_count"

if [ "$initial_fastq_count" -ne "$final_fastq_count" ]; then
    echo "ERROR: FASTQ file count mismatch! Initial: $initial_fastq_count, Final: $final_fastq_count"
    exit 1
fi

echo "Operation completed successfully at $(date)"
echo "Log file: $log_file"
