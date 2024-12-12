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

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
OUTPUT_DIR="${HOME}/data/100303Bel"
FASTQ_DIR="${OUTPUT_DIR}/fastq"
SAMPLES=(
    "WT_G1_37C_Nucleosomes:SRR034477:GSM442537"
    "WT_G2_ORC_rep1:SRR034475:GSM424494"
    "WT_G2_ORC_rep2:SRR034476:GSM424494"
    "orc1-161_G2_37C_Nucleosomes_rep1:SRR034473:GSM424493"
    "orc1-161_G2_37C_Nucleosomes_rep2:SRR034474:GSM424493"
    "WT_G2_37C_Nucleosomes_rep1:SRR034471:GSM424492"
    "WT_G2_37C_Nucleosomes_rep2:SRR034472:GSM424492"
    "WT_Asynchronous_23C_Nucleosomes:SRR034470:GSM424491"
)

# Check for required tools
setup_sratools "prefetch"
setup_sratools "fasterq-dump"

#-------------------------------------------------------------------------------
# Main Process
#-------------------------------------------------------------------------------
# Create output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$FASTQ_DIR"

# Extract accessions from SAMPLES array
ACCESSIONS=()
for sample in "${SAMPLES[@]}"; do
    IFS=':' read -r _ acc _ <<< "$sample"
    ACCESSIONS+=("$acc")
done

# Download and convert each SRA file
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
