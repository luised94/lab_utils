#!/bin/bash
################################################################################
# Clean up intermediates results from pipeline completion tests
################################################################################
# Purpose: Run simple bash commands to remove files generated by exploratory scripts. This was done after initial tests with dataset from V5 set of experiments yielded inconclusive results. The second dataset looked cleaner (especially input samples) by inspecting the bigwig files.
# Usage: I just run them on the command line. No need to turn into super robust script as long as the commands are run from the proper directory. This serves more as documentation for a one-time-operation likely.
# DEPENDENCIES: bash
# OUTPUT: Clean pipeline completion testing
# AUTHOR: LEMR
# DATE: 2025-04-17
################################################################################
# DOCUMENTATION ONLY - DO NOT EXECUTE
# This file contains example code for reference purposes only

echo "ERROR: This file contains documentation code only and is not meant to be executed."
echo "Please open this file in a text editor to view the code examples."
return 2>/dev/null || exit 1

# Check if this script is being executed directly (rather than sourced)
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: This is a documentation-only script and should not be executed."
    echo "It contains example code snippets for reference purposes only."
    echo "Please use 'source $(basename "${0}")' if you need to load the functions."
    exit 1
fi

# Can be adjusted
TARGET_DIRECTORY=~/data/preprocessing_test
cd "${TARGET_DIRECTORY}"  || {
  echo "Directory does not exists. Exiting"
  exit 2
}

# Confirm that all extensions are included by verifying directory before running
# Could have assigned the result of this command to the EXTENSIONS variable but I would like to retain the plots that I have generated so far
find . -type f | awk -F. '{print $NF}' | sort | uniq

EXTENSIONS=(
  "fastq"
  #"sizes"
  "txt"
  "bam"
  "bai"
  "log"
  "r"
  "xls"
  "bed"
  "bw"
  "narrowPeak"
  "broadPeak"
  "gappedPeak"
)

# Run this first to ensure the proper files are found.
for EXTENSION in "${EXTENSIONS[@]}"; do
  formatted_extension="*.${EXTENSION}"
  echo "Finding files for: ${formatted_extension}"
  find ${TARGET_DIRECTORY} -type f -name "${formatted_extension}"
done

for EXTENSION in "${EXTENSIONS[@]}"; do
  formatted_extension="*.${EXTENSION}"
  echo "Finding files for: ${formatted_extension}"
  find ${TARGET_DIRECTORY} -type f -name "${formatted_extension}" -delete
done
