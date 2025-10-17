#!/bin/bash

# Test script for manifest consumption logic
# Simulates what SLURM array tasks would see

# Configuration
EXPERIMENT_DIR="$HOME/data/250930Bel"
MANIFEST_FILE="$EXPERIMENT_DIR/documentation/paired_reads_manifest.tsv"

echo "=== Manifest Consumption Test ==="
echo "Manifest file: $MANIFEST_FILE"
echo ""

# Validate manifest exists
if [[ ! -f "$MANIFEST_FILE" ]]; then
    echo "ERROR: Manifest file not found: $MANIFEST_FILE"
    exit 1
fi

# Count total pairs (excluding header line)
total_pairs=$(tail -n +2 "$MANIFEST_FILE" | wc -l)
echo "Total pairs found: $total_pairs"
echo ""

# Show manifest header
echo "Manifest structure:"
head -n 1 "$MANIFEST_FILE"
echo ""

# Test first 3 task IDs to verify parsing
echo "=== Testing Task ID Parsing ==="
for simulated_task_id in 1 2 3; do
    echo "--- Task ID: $simulated_task_id ---"

    # Validate task ID is within range
    if [[ $simulated_task_id -gt $total_pairs ]]; then
        echo "ERROR: Task ID $simulated_task_id exceeds total pairs ($total_pairs)"
        continue
    fi

    # Extract the line for this task (header is line 1, so add 1 to task ID)
    pair_line=$(sed -n "$((simulated_task_id + 1))p" "$MANIFEST_FILE")

    # Parse TSV columns into variables
    IFS=$'\t' read -r sample_id lane_num read1_path read2_path <<< "$pair_line"

    # Display parsed values
    echo "  Sample ID:  $sample_id"
    echo "  Lane:       $lane_num"
    echo "  Read1:      $read1_path"
    echo "  Read2:      $read2_path"

    # Verify files exist (optional validation)
    if [[ -f "$read1_path" ]]; then
        echo "  û Read1 file exists"
    else
        echo "  ? Read1 file NOT FOUND"
    fi

    if [[ -f "$read2_path" ]]; then
        echo "  û Read2 file exists"
    else
        echo "  ? Read2 file NOT FOUND"
    fi

    echo ""
done

echo "=== Test complete ==="
echo "If parsing looks correct, ready to integrate with SLURM array job"
