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

# Create target directory if it doesn't exist
# Think leave commented since I want the script to tell me if it doesnt exist just in case I forgot to sync or made a mistake.
#mkdir -p "$TARGET_DIR"

# Check if target directory exists and is writable
if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: Target directory '$TARGET_DIR' does not exist"
    exit 1
fi

if [ ! -w "$TARGET_DIR" ]; then
    echo "Error: Target directory '$TARGET_DIR' is not writable"
    exit 1
fi

# Move files in the specified range
for i in $(seq $START_NUM $END_NUM); do
    filename="${PREFIX}${i}${SUFFIX}"
    if [ -f "$SOURCE_DIR/$filename" ]; then
        #mv "$SOURCE_DIR/$filename" "$TARGET_DIR/"
        echo "Moved: $filename"
    else
        echo "File not found: $filename"
    fi
done

echo "File move operation completed"
