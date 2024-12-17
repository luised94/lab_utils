#!/bin/bash

# Configuration
MACS2_ENV="macs2_env"
OUTDIR="macs2_test_results"
GENOME_SIZE="1.2e7"  # S. cerevisiae genome size

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate $MACS2_ENV

# Create output directory
mkdir -p $OUTDIR

## Test Case 1: Sample with Input Control
echo "Running MACS2 with Input Control..."
macs2 callpeak \
    -t path/to/test_sample.bam \
    -c path/to/input_control.bam \
    -n test_with_control \
    -g $GENOME_SIZE \
    --outdir $OUTDIR \
    -q 0.05

## Test Case 2: Reference Sample without Input
echo "Running MACS2 without Input Control..."
macs2 callpeak \
    -t path/to/reference_sample.bam \
    -n test_no_control \
    -g $GENOME_SIZE \
    --outdir $OUTDIR \
    -q 0.05 \
    --nomodel

# Basic Quality Check
echo "Performing basic quality checks..."
wc -l $OUTDIR/test_with_control_peaks.narrowPeak
wc -l $OUTDIR/test_no_control_peaks.narrowPeak
