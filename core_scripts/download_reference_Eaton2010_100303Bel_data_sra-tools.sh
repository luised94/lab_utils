#!/bin/bash
#===============================================================================
# TITLE: Download Eaton 2010 SRA Data
# DESCRIPTION: Downloads and converts SRA files to FASTQ format using SRA-tools
# AUTHOR: [Your Name]
# DATE: 2024-12-12
# VERSION: 1.0.0
#===============================================================================

#-------------------------------------------------------------------------------
# Setup SRA-tools
#-------------------------------------------------------------------------------
setup_sratools() {
    local tool=$1
    
    # First check if tool is already in PATH
    if command -v "$tool" &> /dev/null; then
        return 0
    fi

    # Check if toolkit is installed but not in PATH
    TOOLKIT_DIR=$(ls ${HOME} | grep sratoolkit | grep -v tar.gz)
    if [[ -n "$TOOLKIT_DIR" ]]; then
        echo "Found SRA-tools in ${HOME}/${TOOLKIT_DIR}"
        export PATH="${PATH}:${HOME}/${TOOLKIT_DIR}/bin"
        echo "Added SRA-tools to PATH for current session"
        
        # Verify tool is now available
        if command -v "$tool" &> /dev/null; then
            return 0
        fi
    fi

    # If we get here, tool is not installed
    echo "Error: $tool not found"
    echo "Please run ~/lab_utils/core_scripts/install_ncbi_sra-tools_cli.sh first"
    exit 1
}

# Check for required tools
setup_sratools "prefetch"
setup_sratools "fasterq-dump"

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
OUTPUT_DIR="${HOME}/data/Eaton2010"
ACCESSIONS=(
    "SRR034477"
    "SRR034475"
    "SRR034476"
    "SRR034473"
    "SRR034474"
    "SRR034471"
    "SRR034472"
    "SRR034470"
)

#-------------------------------------------------------------------------------
# Main Process
#-------------------------------------------------------------------------------
# Create output directory
mkdir -p "$OUTPUT_DIR"

# Download each SRA file and convert to FASTQ
for acc in "${ACCESSIONS[@]}"; do
    echo "Processing $acc..."
    
    echo "Downloading $acc..."
    if ! prefetch "$acc" -O "$OUTPUT_DIR"; then
        echo "Error: Failed to download $acc"
        continue
    fi
    
    echo "Converting $acc to FASTQ..."
    if ! fasterq-dump "$OUTPUT_DIR/$acc/$acc.sra" -O "$OUTPUT_DIR" -p; then
        echo "Error: Failed to convert $acc to FASTQ"
        continue
    fi
    
    echo "Successfully processed $acc"
done

echo "All processing completed."
