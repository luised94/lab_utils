#!/bin/bash
#===============================================================================
# TITLE: Download GEO Dataset Files
# DESCRIPTION: Downloads specific files from GEO dataset GSE242131
# DATE: 2024-01-17
# VERSION: 1.0.0
#===============================================================================
set -euo pipefail
IFS=$'\n\t'

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
DOWNLOAD_DIR="${HOME}/downloads/GSE242131"
LOG_FILE="${DOWNLOAD_DIR}/download.log"
TIMESTAMP=$(date '+%Y%m%d_%H%M%S')

# File definitions
declare -A files=(
    ["CMBSs"]="GSE242131_17618_CMBSs_chr_coord_rank_082923_1.txt.gz"
    ["MNase"]="GSE242131_processed_data_for_MNase_seq_size_range_151bpto200bp_083023_1.txt.gz"
    ["script"]="GSE242131_script.tar.gz"
)

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------
log_message() {
    local message="$1"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $message" | tee -a "$LOG_FILE"
}

check_command() {
    if ! command -v "$1" &>/dev/null; then
        log_message "Error: $1 is not installed"
        exit 1
    fi
}

download_file() {
    local url="$1"
    local output="$2"
    local attempt=1
    local max_attempts=3
    
    while [ $attempt -le $max_attempts ]; do
        log_message "Downloading $output (Attempt $attempt/$max_attempts)"
        if wget --no-verbose --show-progress --continue "$url" -O "$output"; then
            log_message "Successfully downloaded $output"
            return 0
        fi
        ((attempt++))
        [ $attempt -le $max_attempts ] && sleep 5
    done
    
    log_message "Failed to download $output after $max_attempts attempts"
    return 1
}

verify_file() {
    local file="$1"
    if [ -f "$file" ]; then
        local size=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null)
        if [ "$size" -gt 0 ]; then
            return 0
        fi
    fi
    return 1
}

#-------------------------------------------------------------------------------
# Main Script
#-------------------------------------------------------------------------------
main() {
    # Check required commands
    check_command wget
    
    # Create download directory
    mkdir -p "$DOWNLOAD_DIR"
    log_message "Created download directory: $DOWNLOAD_DIR"
    
    # Download each file
    local base_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE242nnn/GSE242131/suppl"
    local success=true
    
    for key in "${!files[@]}"; do
        local filename="${files[$key]}"
        local output="${DOWNLOAD_DIR}/${filename}"
        local url="${base_url}/${filename}"
        
        if ! download_file "$url" "$output"; then
            success=false
            continue
        fi
        
        if ! verify_file "$output"; then
            log_message "Error: Failed to verify $filename"
            success=false
            continue
        fi
    done
    
    if [ "$success" = true ]; then
        log_message "All files downloaded successfully"
        return 0
    else
        log_message "Some downloads failed. Check the log file for details"
        return 1
    fi
}

# Execute main function
main

