#!/bin/bash
#===============================================================================
# TITLE: Download and Process Reference Genomes
# DESCRIPTION: Downloads and processes reference genomes using NCBI datasets
# AUTHOR: [Your Name]
# DATE: 2024-12-12
# VERSION: 1.0.0
#===============================================================================

set -euo pipefail

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
DOWNLOAD_DIR="$HOME/data/REFGENS"
LOG_DIR="$HOME/data/REFGENS/logs"
ACCESSIONS=(
    "GCF_000146045.2"  # S. cerevisiae S288C
    "GCF_000001405.40" # Human
    "GCF_000005845.2"  # E. coli
    "GCA_002163515.1"  # S cerevisiae W303
)

#-------------------------------------------------------------------------------
# Setup directories
#-------------------------------------------------------------------------------
mkdir -p "$DOWNLOAD_DIR"
mkdir -p "$LOG_DIR"

#-------------------------------------------------------------------------------
# Download genomes
#-------------------------------------------------------------------------------
echo "Starting genome downloads..."

for accession in "${ACCESSIONS[@]}"; do
    echo "Processing accession: $accession"
    
    # Create timestamped directory
    timestamp=$(date +"%Y%m%d_%H%M%S")
    target_dir="${DOWNLOAD_DIR}/${accession}_${timestamp}"
    mkdir -p "$target_dir"
    
    echo "Downloading genome for $accession to $target_dir"
    if ! datasets download genome accession "$accession" \
        --include genome,rna,cds,protein,gff3,gtf \
        --filename "${target_dir}/${accession}.zip"; then
        echo "Download failed for accession: $accession" >&2
        continue
    fi
    
    echo "Extracting files..."
    unzip "${target_dir}/${accession}.zip" -d "$target_dir"
    rm "${target_dir}/${accession}.zip"
    
    echo "Download and extraction complete for: $accession"
done

#-------------------------------------------------------------------------------
# Build Bowtie2 index (commented out)
#-------------------------------------------------------------------------------
echo "Finding genome files for indexing..."
mapfile -t genome_paths < <(find "${DOWNLOAD_DIR}" -type f -name "GCF*.fna" -o -name "GCA*.fna")

for genome_path in "${genome_paths[@]}"; do
    echo "Would index: $genome_path"
    #bowtie2-build "$genome_path" "${genome_path%.fna}_index"
done

echo "All operations completed successfully."
