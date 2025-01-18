#!/bin/bash
#===============================================================================
# TITLE: Download GEO Dataset Files
# DESCRIPTION: Downloads specific files from GEO dataset GSE242131 with timestamps
# DATE: 2024-01-17
# VERSION: 1.0.0
#===============================================================================
set -euo pipefail
IFS=$'\n\t'

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
FEATURE_DIR="${HOME}/data/feature_files"
DATE_PREFIX=$(date '+%Y%m%d')

# File definitions
declare -A files=(
    ["CMBSs"]="GSE242131_17618_CMBSs_chr_coord_rank_082923_1.txt.gz"
    ["MNase"]="GSE242131_processed_data_for_MNase_seq_size_range_151bpto200bp_083023_1.txt.gz"
    ["script"]="GSE242131_script.tar.gz"
)

#-------------------------------------------------------------------------------
# Main Script
#-------------------------------------------------------------------------------
# Create feature directory if it doesn't exist
mkdir -p "$FEATURE_DIR"
echo "Created directory: $FEATURE_DIR"

# Download each file
base_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242131/suppl"

for key in "${!files[@]}"; do
    filename="${files[$key]}"
    extension="${filename##*.}"
    basename="${filename%.*}"
    # For .tar.gz files, we need to handle double extension
    if [[ "$filename" == *.tar.gz ]]; then
        extension="tar.gz"
        basename="${filename%.tar.gz}"
    fi
    
    timestamped_name="${DATE_PREFIX}_${basename}.${extension}"
    output="${FEATURE_DIR}/${timestamped_name}"
    
    echo "Downloading: ${filename}"
    echo "Saving as: ${timestamped_name}"
    
    if wget --no-verbose --show-progress "$base_url/$filename" -O "$output"; then
        echo "Successfully downloaded: ${timestamped_name}"
    else
        echo "Failed to download: ${filename}"
        exit 1
    fi
done

echo "All downloads completed successfully"
