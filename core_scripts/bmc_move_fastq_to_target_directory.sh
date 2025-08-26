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

# Configuration variables: Adjust all of these as necessary
#SOURCE_DIR="${1-.}"
SOURCE_DIR="$HOME/data/250715Bel/"
TARGET_DIR="$HOME/data/250714Bel/fastq"
START_NUM=219001
END_NUM=219030

# Expand and display directory paths
SOURCE_DIR=$(readlink -f "$SOURCE_DIR")
TARGET_DIR=$(readlink -f "$TARGET_DIR")
echo "Source directory: $SOURCE_DIR"
echo "Target directory: $TARGET_DIR"
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
missing_id=()
all_files=()

for id in $(seq "$START_NUM" "$END_NUM"); do
  files_for_id=()
  # Null-delimited from find, handles spaces/newlines
  while IFS= read -r -d $'\0' f; do
    files_for_id+=("$f")
  done < <(find "$SOURCE_DIR" -type f -name "*${id}*.fastq" -print0)

  if (( ${#files_for_id[@]} == 0 )); then
    missing_id+=("$id")
  else
    printf 'Found: %s (%d files)\n' "$id" "${#files_for_id[@]}"
    all_files+=("${files_for_id[@]}")
  fi
done

if (( ${#missing_id[@]} > 0 )); then
  printf 'Missing IDs (%d): %s\n' "${#missing_id[@]}" "${missing_id[*]}"
  exit 1
fi

# Move only if none missing
if (( ${#all_files[@]} )); then
  # dry run:
  printf 'Would move: %s\n' "${all_files[@]}"
  #echo "mv -t "$TARGET_DIR" -- "${all_files[@]}""
  #mv -t "$TARGET_DIR" -- "${all_files[@]}"
  printf 'Moved total: %d files\n' "${#all_files[@]}"
else
  echo 'Nothing to move.'
fi
