#!/bin/bash

set -euo pipefail

# Script Name: refactored_genome_processing.sh
# Description: Download, reorganize, and index reference genomes using NCBI datasets
# Usage: ./refactored_genome_processing.sh [OPTIONS]

# Global variables
DOWNLOAD_DIR="$HOME/data/REFGENS"
LOG_DIR="$HOME/data/REFGENS/logs"

# Function to download genomes using datasets
download_genomes() {
    local accessions=("$@")
    mkdir -p "$DOWNLOAD_DIR"

    for accession in "${accessions[@]}"; do
        echo "Downloading genome for accession: $accession"
        local timestamp=$(date +"%Y%m%d_%H%M%S")
        local target_dir="${DOWNLOAD_DIR}/${accession}_${timestamp}"
        mkdir -p "$target_dir"
        if ! datasets download genome accession "$accession" --include genome,rna,cds,protein,gff3,gtf --filename "${target_dir}/${accession}.zip"; then
            echo "Download failed for accession: $accession. Check network or accession validity." >&2
            continue
        fi

        unzip "${target_dir}/${accession}.zip" -d "$target_dir"
        rm "${target_dir}/${accession}.zip"

        echo "Download complete for: $accession"
        echo "Downloaded to: $target_dir"
    done
}

# Function to build Bowtie2 index
build_bowtie2_index() {
    local refgenome_dir="$HOME/data/REFGENS"
    local genome_paths
    mapfile -t genome_paths < <(find "${refgenome_dir}" -type f -name "GCF*.fna" -o -name "GCA*.fna")

    for genome_path in "${genome_paths[@]}"; do
        echo "Starting indexing for $genome_path"
        #bowtie2-build "$genome_path" "${genome_path%.fna}_index"
        echo "Indexing completed for $genome_path"
    done
}

# Main function
main() {
    local accessions=(
        "GCF_000146045.2"  # S. cerevisiae S288C
        "GCF_000001405.40" # Human
        "GCF_000005845.2"  # E. coli
        "GCA_002163515.1"  # S cerevisaie W303
    )
    download_genomes "${accessions[@]}"
    #build_bowtie2_index
    echo "All operations completed successfully."
}

# Run the main function
main "$@"
