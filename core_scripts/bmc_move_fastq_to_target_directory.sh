#!/bin/bash
################################################################################
# Move from particular directory to target directory
# Author: Luis | Date: 2025-08-25 | Version: 1.0.0
################################################################################
# PURPOSE:
#   Move files to fastq directory to group with respective metadata
# USAGE:
#   Run from command line line by line.
# DEPENDENCIES:
#   Data from bmc
# OUTPUTS:
#   Moves fastq files that follow prefix-num-suffix pattern to target directory
################################################################################

# Configuration variables
SOURCE_DIR="${1-.}"
# Adjust all of these as necessary
TARGET_DIR="/path/to/target/directory"
PREFIX="250715Bel_D25-"
SUFFIX="-1_NA_sequence.fastq"
START_NUM=219001
END_NUM=219037

# Expand and display directory paths
SOURCE_DIR=$(realpath "$SOURCE_DIR")
TARGET_DIR=$(realpath "$TARGET_DIR")
echo "Source directory: $SOURCE_DIR"
echo "Target directory: $TARGET_DIR"
echo "File pattern: ${PREFIX}*${SUFFIX}"
echo "Number range: $START_NUM to $END_NUM"
echo

# Check if source directory exists and is readable
if [ ! -d "$SOURCE_DIR" ]; then
    echo "Error: Source directory '$SOURCE_DIR' does not exist"
    exit 1
fi

if [ ! -r "$SOURCE_DIR" ]; then
    echo "Error: Source directory '$SOURCE_DIR' is not readable"
    exit 1
fi

# Check if target directory exists and is writable
if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: Target directory '$TARGET_DIR' does not exist"
    exit 1
fi

if [ ! -w "$TARGET_DIR" ]; then
    echo "Error: Target directory '$TARGET_DIR' is not writable"
    exit 1
fi

# First pass: Check that all files exist
echo "Checking file availability..."
missing_files=()
for i in $(seq $START_NUM $END_NUM); do
    filename="${PREFIX}${i}${SUFFIX}"
    if [ ! -f "$SOURCE_DIR/$filename" ]; then
        missing_files+=("$filename")
        echo "Missing: $filename"
    else
        echo "Found: $filename"
    fi
done

# Exit if any files are missing
if [ ${#missing_files[@]} -gt 0 ]; then
    echo
    echo "Error: ${#missing_files[@]} files are missing. Aborting operation."
    echo "Verify START_NUM, END_NUM, SOURCE_DIR, PREFIX, SUFFIX."
    exit 1
fi

echo
echo "All files found. Proceeding with move operation..."
echo

# Second pass: Move all files
#for i in $(seq $START_NUM $END_NUM); do
#    filename="${PREFIX}${i}${SUFFIX}"
#    mv "$SOURCE_DIR/$filename" "$TARGET_DIR/"
#    echo "Moved: $filename"
#done

echo
echo "File move operation completed successfully"
