#!/bin/bash
#===============================================================================
# TITLE: Download Eaton 2010 SRA Data
# DESCRIPTION: Downloads and converts SRA files to FASTQ format using SRA-tools
# AUTHOR: [Your Name]
# DATE: 2024-12-12
# VERSION: 1.0.0
#===============================================================================
<<<<<<< HEAD
=======

>>>>>>> 48e7cdd (- Add Eaton2010 data download script using curl and sra-tools.)
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

# Download each SRA file and convert to FASTQ
>>>>>>> 48e7cdd (- Add Eaton2010 data download script using curl and sra-tools.)
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

#-------------------------------------------------------------------------------
# Cleanup and Consolidation
#-------------------------------------------------------------------------------
echo "Organizing files..."

# Move all FASTQ files to FASTQ directory
find "$OUTPUT_DIR" -type f -name "*.fastq" -exec mv {} "$FASTQ_DIR" \;

# Remove all non-FASTQ files
find "$OUTPUT_DIR" -type f ! -name "*.fastq" -delete

# Remove empty directories
find "$OUTPUT_DIR" -type d -empty -delete

echo "Processing completed. FASTQ files are in ${FASTQ_DIR}"

echo "Starting consolidation process..."
ls -l "$FASTQ_DIR"

read -p "Continue with these files? (y/n): " confirm
if [[ ! $confirm =~ ^[Yy]$ ]]; then
    echo "Operation cancelled"
    exit 0
fi
#
# Group samples by GSM ID
declare -A gsm_groups
for sample in "${SAMPLES[@]}"; do
    IFS=':' read -r name srr gsm <<< "$sample"
    if [[ ! ${gsm_groups[$gsm]} ]]; then
        gsm_groups[$gsm]=""
    fi
    gsm_groups[$gsm]+="${FASTQ_DIR}/${srr}.fastq "
done

# Show consolidation plan
echo -e "\nPlanned consolidation:"
for gsm in "${!gsm_groups[@]}"; do
    echo "GSM: $gsm"
    echo "Files: ${gsm_groups[$gsm]}"
    echo "---"
done

read -p "Proceed with consolidation? (y/n): " proceed
if [[ ! $proceed =~ ^[Yy]$ ]]; then
    echo "Consolidation cancelled"
    exit 0
fi
# Process each GSM group
for gsm in "${!gsm_groups[@]}"; do
    echo "Processing GSM: $gsm"
    
    # Get first SRR for this group for naming
    first_srr=$(echo "${SAMPLES[@]}" | tr ' ' '\n' | grep "$gsm" | head -n1 | cut -d: -f2)
    id=${first_srr:3:6}
    
    files=${gsm_groups[$gsm]}
    output_file="${FASTQ_DIR}/consolidated_${id}_sequence.fastq"
    
    echo "Consolidating files into $output_file..."
    echo "Input files: $files"
    
    if cat $files > "$output_file"; then
        if [ -s "$output_file" ]; then
            original_size=$(du -b $files | awk '{sum += $1} END {print sum}')
            new_size=$(du -b "$output_file" | awk '{print $1}')
            
            echo "Original files total size: $original_size bytes"
            echo "New file size: $new_size bytes"
            
            if [ "$new_size" -gt 0 ]; then
                echo "Removing original files..."
                rm -f $files
                echo "Original files removed"
            else
                echo "Error: Consolidated file is empty"
                rm -f "$output_file"
                continue
            fi
        fi
    else
        echo "Error during consolidation"
        rm -f "$output_file"
        continue
    fi
done

echo "Consolidation completed"
echo -e "\nTo sync data manually, run the following command:\n"
echo "rsync -avzP LOCAL_DIR/ REMOTE_HOST:REMOTE_DIR/"
echo -e "\nReplace:\nLOCAL_DIR with ${HOME}/data/100303Bel\nREMOTE_HOST with your Luria username@luria.mit.edu\nREMOTE_DIR with ~/data/100303Bel"
echo "All consolidation completed successfully"
