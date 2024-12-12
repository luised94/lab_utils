#!/bin/bash
#===============================================================================
# TITLE: Download and Process Reference Genomes
# DESCRIPTION: Downloads and processes reference genomes using NCBI datasets
# AUTHOR: [Your Name]
# DATE: 2024-12-12
# VERSION: 2.0.0
#===============================================================================

set -euo pipefail

#-------------------------------------------------------------------------------
# Setup NCBI datasets
#-------------------------------------------------------------------------------
setup_ncbi_datasets() {
    local tool="datasets"
    
    if command -v "$tool" &> /dev/null; then
        return 0
    fi

    DATASETS_DIR=$(ls "${HOME}" | grep "ncbi-datasets-cli" | grep -v tar.gz)
    if [[ -n "$DATASETS_DIR" ]]; then
        echo "Found NCBI datasets CLI in ${HOME}/${DATASETS_DIR}"
        export PATH="${PATH}:${HOME}/${DATASETS_DIR}/bin"
        
        if command -v "$tool" &> /dev/null; then
            return 0
        fi
    fi

    echo "Error: $tool not found"
    echo "Please run ~/lab_utils/core_scripts/install_ncbi_datasets_cli.sh first"
    exit 1
}

setup_ncbi_datasets
#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
DOWNLOAD_DIR="$HOME/data/REFGENS"
LOG_DIR="$HOME/data/REFGENS/logs"
DOWNLOAD_LOG="${LOG_DIR}/reference_genome_metadata.txt"

# File naming patterns
ORIGINAL_CDS="cds_from_genomic.fna"
CDS_FILE="cds.fna"
GENOMIC_PATTERN="*.fna"
REFGENOME_SUFFIX="_refgenome.fna"
NCBI_DATA="ncbi_dataset"

mkdir -p "$LOG_DIR"

echo -e "\nConfiguration:"
echo "Download directory: $DOWNLOAD_DIR"
echo "Log directory: $LOG_DIR"
echo "Download log file: $DOWNLOAD_LOG"
echo -e "\nAccessions to process:"
printf '%s\n' "${ACCESSIONS[@]}"

read -p "Continue with these settings? (y/n): " confirm
if [[ ! $confirm =~ ^[Yy]$ ]]; then
    echo "Operation cancelled"
    exit 0
fi
#-------------------------------------------------------------------------------
# Download and Process
#-------------------------------------------------------------------------------
ACCESSIONS=(
    "GCF_000146045.2"  # S. cerevisiae S288C
    "GCF_000001405.40" # Human
    "GCF_000005845.2"  # E. coli
    "GCA_002163515.1"  # S cerevisiae W303
)

echo "Starting genome downloads..."

for accession in "${ACCESSIONS[@]}"; do
    echo "Processing accession: $accession"
    target_dir="${DOWNLOAD_DIR}/${accession}"

    echo -e "\nProcessing details for $accession:"
    echo "Target directory: $target_dir"
    echo "Expected files: genome, RNA, CDS, protein, GFF3, GTF"
    
    read -p "Process this accession? (y/n): " process_confirm
    if [[ ! $process_confirm =~ ^[Yy]$ ]]; then
        echo "Skipping $accession"
        continue
    fi

    mkdir -p "$target_dir"
    
    # Log download timestamp
    echo "$(date '+%Y-%m-%d %H:%M:%S') - Downloading $accession" >> "$DOWNLOAD_LOG"
    
    # Download genome
    if ! datasets download genome accession "$accession" \
        --include genome,rna,cds,protein,gff3,gtf \
        --filename "${target_dir}/${accession}.zip"; then
        echo "Download failed for accession: $accession" >&2
        continue
    fi
    
    # Extract files
    unzip "${target_dir}/${accession}.zip" -d "$target_dir"
    rm "${target_dir}/${accession}.zip"
    
    # Extract organism name from metadata
    report_file="${target_dir}/${NCBI_DATA}/data/assembly_data_report.jsonl"
    organism_name=$(grep -m 1 -o '"organismName":"[^"]*' "$report_file" | 
                   cut -d '"' -f4 | 
                   sed 's/ //g')
    
    echo -e "\nExtracted information:"
    echo "Organism name: $organism_name"
    echo "Source directory: $target_dir"
    echo "Target CDS file: ${target_dir}/${CDS_FILE}"
    echo "Target genome file: ${target_dir}/${organism_name}${REFGENOME_SUFFIX}"
    
    read -p "Continue with file reorganization? (y/n): " reorg_confirm
    if [[ ! $reorg_confirm =~ ^[Yy]$ ]]; then
        echo "Skipping reorganization for $accession"
        continue
    fi

    if [ -z "$organism_name" ]; then
        echo "Failed to extract organism name for $accession"
        continue
    fi
    
    echo "Reorganizing files for: $organism_name"
    
    # Move all files to root directory
    find "$target_dir/${NCBI_DATA}/data/" -type f -exec mv -v {} "$target_dir/" \;
    
    # Rename CDS file if exists
    if [ -f "${target_dir}/${ORIGINAL_CDS}" ]; then
        mv -v "${target_dir}/${ORIGINAL_CDS}" "${target_dir}/${CDS_FILE}"
    fi
    
    # Rename genomic file
    find "$target_dir" -name "${GENOMIC_PATTERN}" ! -name "${CDS_FILE}" \
        -exec mv -v {} "${target_dir}/${organism_name}${REFGENOME_SUFFIX}" \;
    
    # Cleanup
    rm -rf "${target_dir}/${NCBI_DATA}"
    
    # Rename directory to organism name
    new_dir="${DOWNLOAD_DIR}/${organism_name}"
    if [[ -d "$new_dir" ]]; then
        echo "Warning: Directory $new_dir already exists"
        read -p "Overwrite existing directory? (y/n): " overwrite_confirm
        if [[ ! $overwrite_confirm =~ ^[Yy]$ ]]; then
            echo "Skipping $accession"
            continue
        fi
        rm -rf "$new_dir"
    fi

    echo "Renaming directory from ${target_dir} to ${new_dir}"
    mv "$target_dir" "$new_dir"
    target_dir="$new_dir"

    # Log metadata
    {
        echo "----------------------------------------"
        echo "Organism: $organism_name"
        echo "Accession: $accession"
        echo "Download Date: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "Download Directory: $new_dir"
        echo "----------------------------------------"
    } >> "${new_dir}/reference_genome_metadata.txt"
    echo "Completed processing: $accession"

done

echo "All operations completed successfully."

#-------------------------------------------------------------------------------
# Build Bowtie2 index (commented out)
#-------------------------------------------------------------------------------
#echo "Finding genome files for indexing..."
#mapfile -t genome_paths < <(find "${DOWNLOAD_DIR}" -type f -name "GCF*.fna" -o -name "GCA*.fna")
#
#for genome_path in "${genome_paths[@]}"; do
#    echo "Would index: $genome_path"
#    #bowtie2-build "$genome_path" "${genome_path%.fna}_index"
#done
#
#echo "All operations completed successfully."
