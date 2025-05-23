#!/bin/bash

# Script to download yeast IDR (Intrinsically Disordered Region) FASTA files
# Author: Data Analysis Script
# Date: $(date +%Y-%m-%d)
# Purpose: Reproducible download of yeast IDR datasets for analysis

set -e  # Exit on any error

# Define URLs and output filenames
IUPRED_URL="https://owncloud.mpi-cbg.de/index.php/s/erQbBdBPUfLJbBA/download?path=%2Fiupred_based&files=yeast_IDRs.fa"
ALPHAFOLD_URL="https://owncloud.mpi-cbg.de/index.php/s/erQbBdBPUfLJbBA/download?path=%2Falphafold2_based&files=UP000002311_559292_YEAST_v4_min10_IDRs.fa"

IUPRED_OUTPUT="yeast_IDRs_iupred.fa"
ALPHAFOLD_OUTPUT="yeast_IDRs_alphafold2.fa"

# Create data directory if it doesn't exist
DATA_DIR="data"
mkdir -p "$DATA_DIR"

echo "=== Downloading Yeast IDR FASTA Files ==="
echo "Data will be saved to: $DATA_DIR/"
echo

# Function to download with error handling
download_file() {
    local url=$1
    local output=$2
    local description=$3
    
    echo "Downloading $description..."
    echo "URL: $url"
    echo "Output: $DATA_DIR/$output"
    
    if curl -L -o "$DATA_DIR/$output" "$url"; then
        echo "û Successfully downloaded $output"
        
        # Display file info
        if [[ -f "$DATA_DIR/$output" ]]; then
            file_size=$(du -h "$DATA_DIR/$output" | cut -f1)
            echo "  File size: $file_size"
            
            # Count sequences if it's a FASTA file
            if command -v grep >/dev/null 2>&1; then
                seq_count=$(grep -c "^>" "$DATA_DIR/$output" 2>/dev/null || echo "Unable to count")
                echo "  Number of sequences: $seq_count"
            fi
        fi
    else
        echo "? Failed to download $output"
        return 1
    fi
    echo
}

# Download IUPred-based IDRs
download_file "$IUPRED_URL" "$IUPRED_OUTPUT" "IUPred-based yeast IDRs"

# Download AlphaFold2-based IDRs
download_file "$ALPHAFOLD_URL" "$ALPHAFOLD_OUTPUT" "AlphaFold2-based yeast IDRs"

echo "=== Download Summary ==="
echo "Files downloaded to $DATA_DIR/:"
ls -lh "$DATA_DIR"/*.fa 2>/dev/null || echo "No FASTA files found"

echo
echo "=== Usage Notes ==="
echo "- IUPred-based file: Contains IDRs predicted using IUPred method"
echo "- AlphaFold2-based file: Contains IDRs identified from AlphaFold2 structures"
echo "- Both files are in FASTA format and ready for downstream analysis"

echo
echo "Script completed successfully!"
