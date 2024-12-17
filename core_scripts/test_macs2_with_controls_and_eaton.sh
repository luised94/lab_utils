#!/bin/bash

# Exit on error
set -e

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

# Function to handle errors
error_exit() {
    echo -e "${RED}ERROR: $1${NC}" >&2
    exit 1
}

# Function to validate file existence
validate_file() {
    if [ ! -f "$1" ]; then
        error_exit "File not found: $1"
    fi
}

# Configuration
MACS2_ENV="macs2_env"
OUTDIR="macs2_test_results"
GENOME_SIZE="1.2e7"  # S. cerevisiae genome size
TEST_SAMPLE="path/to/test_sample.bam"
INPUT_CONTROL="path/to/input_control.bam"
REF_SAMPLE="path/to/reference_sample.bam"

# Source conda/MACS2 setup
if ! source ~/conda_macs2_setup.sh; then
    error_exit "Failed to setup MACS2 environment"
fi

# Verify MACS2 functionality
if ! macs2 callpeak -h &>/dev/null; then
    error_exit "MACS2 callpeak command not functioning"
fi

# Validate input files
validate_file "$TEST_SAMPLE"
validate_file "$INPUT_CONTROL"
validate_file "$REF_SAMPLE"

# Create output directory
mkdir -p "$OUTDIR" || error_exit "Failed to create output directory"

# Function to check output file creation
check_output() {
    if [ ! -f "$1" ]; then
        error_exit "Output file not created: $1"
    fi
    echo -e "${GREEN}Successfully created: $1${NC}"
}

## Test Case 1: Sample with Input Control
echo "Running MACS2 with Input Control..."
macs2 callpeak \
    -t "$TEST_SAMPLE" \
    -c "$INPUT_CONTROL" \
    -n test_with_control \
    -g "$GENOME_SIZE" \
    --outdir "$OUTDIR" \
    -q 0.05 || error_exit "MACS2 failed on test with control"

# Verify output
check_output "$OUTDIR/test_with_control_peaks.narrowPeak"

## Test Case 2: Reference Sample without Input
echo "Running MACS2 without Input Control..."
macs2 callpeak \
    -t "$REF_SAMPLE" \
    -n test_no_control \
    -g "$GENOME_SIZE" \
    --outdir "$OUTDIR" \
    -q 0.05 \
    --nomodel || error_exit "MACS2 failed on test without control"

# Verify output
check_output "$OUTDIR/test_no_control_peaks.narrowPeak"

# Basic Quality Check
echo "Performing basic quality checks..."
for peak_file in "$OUTDIR"/*.narrowPeak; do
    echo "Checking $peak_file:"
    wc -l "$peak_file"
    if [ ! -s "$peak_file" ]; then
        error_exit "Peak file is empty: $peak_file"
    fi
done

echo -e "${GREEN}All tests completed successfully${NC}"
