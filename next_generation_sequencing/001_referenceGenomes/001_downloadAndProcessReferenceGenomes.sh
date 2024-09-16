#!/bin/bash

set -euo pipefail

# Script Name: refactored_genome_processing.sh
# Description: Download, reorganize, and index reference genomes from FTP
# Usage: ./refactored_genome_processing.sh [OPTIONS]

# Global variables
DOWNLOAD_DIR="$HOME/data/REFGENS"
LOG_DIR="$HOME/data/REFGENS/logs"

# Function to download genomes from FTP
download_genomes() {
    local urls=("$@")
    mkdir -p "$DOWNLOAD_DIR"

    for url in "${urls[@]}"; do
        echo "Downloading genome from: $url"
        local organism_name=$(echo "$url" | awk -F'/' '{print $(NF-1)}')
        local timestamp=$(date +"%Y%m%d_%H%M%S")
        local target_dir="${DOWNLOAD_DIR}/${organism_name}_${timestamp}"
        #echo "wget --recursive --no-parent --no-host-directories --cut-dirs=7 --directory-prefix=$target_dir $url"
        if ! wget --recursive --no-parent --no-host-directories --cut-dirs=7 --directory-prefix="$target_dir" "$url"; then
            echo "Download failed for URL: $url. Check network or URL validity." >&2
            continue
        fi

        echo "Download complete for: $organism_name"
        echo "Downloaded to: $target_dir"
    done
}

# Function to reorganize genome directories
#reorganize_genome_dirs() {
#    local dirs=("$@")
#    for dir in "${dirs[@]}"; do
#        echo "Reorganizing $dir"
#        local organism_name=$(basename "$dir" | cut -d'_' -f1)
#        echo "Organism name is $organism_name"
#
#        # Move all files to the root of the directory
#        find "$dir" -type f -exec mv -v {} "$dir/" \;
#
#        # Rename files as needed
#        if [[ -f "$dir/cds_from_genomic.fna" ]]; then
#            mv -v "$dir/cds_from_genomic.fna" "$dir/cds.fna"
#        fi
#        
#        local genomic_file=$(find "$dir" -name "*_genomic.fna" | head -n 1)
#        if [[ -n "$genomic_file" ]]; then
#            mv -v "$genomic_file" "$dir/${organism_name}_refgenome.fna"
#        fi
#
#        echo "Reorganization completed for $organism_name"
#    done
#}
#
## Function to rename S288C genome headers (if needed)
#rename_s288c_headers() {
#    local refgenome_dir="$HOME/data/REFGENS"
#    local path_to_s288c_genome
#    path_to_s288c_genome=$(find "$refgenome_dir" -type f -name "*S288C*_refgenome.fna")
#    
#    if [[ -n "$path_to_s288c_genome" ]]; then
#        local backup_file="${path_to_s288c_genome%_refgenome.fna}_backup.fna"
#
#        cp "$path_to_s288c_genome" "$backup_file"
#        awk '/^>/ {gsub(/chromosome/, "chr", $6); printf(">%s%s\n", $6, $7)} !/^>/ {print $0}' "$backup_file" | sed 's/,.*//' > "$path_to_s288c_genome"
#        echo "S288C genome headers renamed"
#    else
#        echo "S288C genome not found, skipping header renaming"
#    fi
#}
#
## Function to build Bowtie2 index
#build_bowtie2_index() {
#    local refgenome_dir="$HOME/data/REFGENS"
#    local genome_paths
#    mapfile -t genome_paths < <(find "${refgenome_dir}" -type f -name "*_refgenome.fna")
#
#    for genome_path in "${genome_paths[@]}"; do
#        echo "Starting indexing for $genome_path"
#        bowtie2-build "$genome_path" "${genome_path%_refgenome.fna}_index"
#        echo "Indexing completed for $genome_path"
#    done
#}

# Main function
main() {
    local urls=(
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/"
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/"
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/"
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/163/515/GCA_002163515.1_ASM216351v1/"
    )
    download_genomes "${urls[@]}"
#    local download_dirs
#    mapfile -t download_dirs < <(find "$DOWNLOAD_DIR" -maxdepth 1 -type d)
#    reorganize_genome_dirs "${download_dirs[@]}"
#    rename_s288c_headers
#    build_bowtie2_index
#    echo "All operations completed successfully."
}

# Run the main function
main "$@"
