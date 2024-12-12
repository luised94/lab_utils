#!/bin/bash
#===============================================================================
# TITLE: Download Feature Data
# DESCRIPTION: Downloads feature data from Rossi 2021 paper and additional sources
# AUTHOR: [Your Name]
# DATE: 2024-12-12
# VERSION: 1.0.0
#===============================================================================

set -euo pipefail

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
FEATURE_DIR="${HOME}/data/feature_files"
DOWNLOAD_LOG="${FEATURE_DIR}/feature_download_metadata.txt"

# Repository configuration
ROSSI_REPO_URL="https://github.com/CEGRcode/2021-Rossi_Nature.git"
ROSSI_BRANCH="main"
ROSSI_DEPTH=1
ROSSI_DIR="rossi_2021"

# Additional data URLs
HAWKINS_TIMING_URL="https://ars.els-cdn.com/content/image/1-s2.0-S2211124713005834-mmc2.xlsx"
EATON_ORC_BED_URL="https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM424nnn/GSM424494/suppl/GSM424494_wt_G2_orc_chip_combined.bed.gz"

#-------------------------------------------------------------------------------
# Git validation
#-------------------------------------------------------------------------------
if ! command -v git &>/dev/null; then
    echo "Error: git is not installed"
    exit 1
fi

#-------------------------------------------------------------------------------
# Directory setup and validation
#-------------------------------------------------------------------------------
mkdir -p "$FEATURE_DIR"

echo -e "\nDirectory structure:"
echo "Feature directory: $FEATURE_DIR"

read -p "Continue with these directories? (y/n): " confirm
if [[ ! $confirm =~ ^[Yy]$ ]]; then
    echo "Operation cancelled"
    exit 0
fi

#-------------------------------------------------------------------------------
# Download Rossi data
#-------------------------------------------------------------------------------
target_dir="${FEATURE_DIR}/${ROSSI_DIR}"
echo "Downloading Rossi 2021 data to: $target_dir"

if [[ -d "$target_dir" ]]; then
    echo "Warning: Directory already exists: $target_dir"
    read -p "Overwrite? (y/n): " overwrite
    if [[ ! $overwrite =~ ^[Yy]$ ]]; then
        exit 0
    fi
    rm -rf "$target_dir"
fi

if ! git clone --depth="$ROSSI_DEPTH" -b "$ROSSI_BRANCH" "$ROSSI_REPO_URL" "$target_dir"; then
    echo "Error: Failed to clone repository"
    exit 1
fi

# Log Rossi download
{
    echo "----------------------------------------"
    echo "Repository: Rossi 2021"
    echo "URL: $ROSSI_REPO_URL"
    echo "Branch: $ROSSI_BRANCH"
    echo "Download Date: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "----------------------------------------"
} >> "$DOWNLOAD_LOG"

#-------------------------------------------------------------------------------
# Download additional data
#-------------------------------------------------------------------------------
feature_folder="${FEATURE_DIR}/${ROSSI_DIR}/02_References_and_Features_Files"

# Download Hawkins timing data
echo "Downloading Hawkins timing data..."
if ! curl -L -o "${feature_folder}/hawkins-origins-timing.xlsx" "$HAWKINS_TIMING_URL"; then
    echo "Error: Failed to download Hawkins timing data"
    exit 1
fi

# Download Eaton ORC bed data
echo "Downloading Eaton ORC bed data..."
if ! curl -L -o "${feature_folder}/G2_orc_chip.bed.gz" "$EATON_ORC_BED_URL"; then
    echo "Error: Failed to download Eaton ORC bed data"
    exit 1
fi

# Extract bed file
gunzip -f "${feature_folder}/G2_orc_chip.bed.gz"

# Log additional downloads
{
    echo "----------------------------------------"
    echo "Additional Data Downloads:"
    echo "Hawkins timing data: $HAWKINS_TIMING_URL"
    echo "Eaton ORC bed data: $EATON_ORC_BED_URL"
    echo "Download Date: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "----------------------------------------"
} >> "$DOWNLOAD_LOG"

echo "All downloads completed successfully"
