#!/bin/bash

# ============================================
# Configuration
# ============================================
FASTQ_DIRECTORY="$HOME/data/250930Bel/fastq"
EXPECTED_LANES_PER_SAMPLE=3
EXPECTED_READ_PAIRS_PER_LANE=2  # R1 and R2
EXPECTED_FILES_PER_SAMPLE=$((EXPECTED_LANES_PER_SAMPLE * EXPECTED_READ_PAIRS_PER_LANE))

# Filename parsing indices (split on _ and -)
# Example: 250930Bel_D25-12496-2_1_sequence.fastq
# Parts: [250930Bel, D25, 12496, 2, 1, sequence.fastq]
SAMPLE_ID_START_IDX=1   # "D25"
SAMPLE_ID_END_IDX=2     # "12496"
LANE_IDX=3              # "2"
READ_PAIR_IDX=4         # "1" or "2"

echo "Looking for fastq files in: $FASTQ_DIRECTORY"
echo "Expected files per sample: $EXPECTED_FILES_PER_SAMPLE"
find "$FASTQ_DIRECTORY" -type f -name "*.fastq" | wc -l
