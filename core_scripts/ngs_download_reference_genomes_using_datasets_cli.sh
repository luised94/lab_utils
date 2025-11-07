#!/bin/bash
set -euo pipefail

# Show usage information
if [[ "${1:-}" == "-h" ]] || [[ "${1:-}" == "--help" ]]; then
    cat << 'EOF'
USAGE: ./ngs_download_reference_genomes_using_datasets_cli.sh

DESCRIPTION:
    Downloads reference genomes from NCBI and organizes them into a flat
    directory structure with standardized naming.

EXAMPLES:
    ./download_reference_genomes.sh

REQUIREMENTS:
    - NCBI datasets CLI tool (datasets command)
    - Sufficient disk space (estimate 10GB per large genome)
    - Internet connection

OUTPUT:
    Files organized in: ~/data/reference_genomes/
    Naming: Organism_Strain_AccessionSanitized_filetype.ext

EOF
    exit 0
fi


#===============================================================================
# Download and Process Reference Genomes
# Downloads genomes from NCBI and organizes into flat directory structure
#===============================================================================

# Configuration - all paths and settings
home_directory="$HOME"
download_base_directory="${home_directory}/data/reference_genomes"
temp_download_directory="${download_base_directory}/temp_downloads"
log_file_path="${download_base_directory}/download_log.txt"

# Accessions to download
accession_list=(
    "GCF_000146045.2"  # S. cerevisiae S288C
    "GCF_000001405.40" # Human
    "GCF_000005845.2"  # E. coli
    "GCA_002163515.1"  # S. cerevisiae W303
)

# File types we expect from NCBI
expected_genome_file="genomic.fna"
expected_cds_file="cds_from_genomic.fna"
expected_protein_file="protein.faa"
expected_gff_file="genomic.gff"
expected_assembly_metadata="assembly_data_report.jsonl"
expected_sequence_metadata="sequence_report.jsonl"
expected_dataset_catalog="dataset_catalog.json"

# Create directories if they don't exist
mkdir -p "$download_base_directory"
mkdir -p "$temp_download_directory"

# Check if datasets tool is available
if ! command -v datasets &> /dev/null; then
    echo "ERROR: 'datasets' command not found"
    echo "Please install NCBI datasets CLI first"
    echo "Installation script: ~/lab_utils/core_scripts/install_ncbi_datasets_cli.sh"
    exit 1
fi

# Show configuration summary
echo "========================================"
echo "Reference Genome Download Configuration"
echo "========================================"
echo "Base directory: $download_base_directory"
echo "Temporary downloads: $temp_download_directory"
echo "Log file: $log_file_path"
echo ""
echo "Accessions to process:"
for accession in "${accession_list[@]}"; do
    echo "  - $accession"
done
echo ""
echo "Total genomes to download: ${#accession_list[@]}"
echo "========================================"
echo ""

# Get user confirmation
read -p "Continue with download? (y/n): " user_confirmation
if [[ ! $user_confirmation =~ ^[Yy]$ ]]; then
    echo "Download cancelled by user"
    exit 0
fi

echo ""
echo "Starting downloads..."
echo "$(date '+%Y-%m-%d %H:%M:%S') - Download started" >> "$log_file_path"

# Process each accession
for accession in "${accession_list[@]}"; do
    echo ""
    echo "========================================"
    echo "Processing: $accession"
    echo "========================================"

    # Create temporary directory for this accession
    accession_temp_directory="${temp_download_directory}/${accession}"
    mkdir -p "$accession_temp_directory"

    # Set download file path
    download_zip_file="${accession_temp_directory}/${accession}.zip"

    echo "Downloading genome data..."
    echo "  Target: $download_zip_file"

    # Download genome with all required data types
    if ! datasets download genome accession "$accession" \
        --include genome,rna,cds,protein,gff3,gtf \
        --filename "$download_zip_file"; then
        echo "ERROR: Download failed for $accession"
        echo "$(date '+%Y-%m-%d %H:%M:%S') - FAILED: $accession" >> "$log_file_path"
        echo "Skipping to next accession..."
        continue
    fi

    echo "Download completed successfully"

    # Verify zip file exists and has content
    if [[ ! -f "$download_zip_file" ]]; then
        echo "ERROR: Download file not found: $download_zip_file"
        echo "$(date '+%Y-%m-%d %H:%M:%S') - FAILED (no file): $accession" >> "$log_file_path"
        continue
    fi

    download_file_size_bytes=$(stat -f%z "$download_zip_file" 2>/dev/null || stat -c%s "$download_zip_file" 2>/dev/null)
    echo "Downloaded file size: $download_file_size_bytes bytes"

    if [[ $download_file_size_bytes -lt 1000 ]]; then
        echo "ERROR: Download file too small, likely failed"
        echo "$(date '+%Y-%m-%d %H:%M:%S') - FAILED (file too small): $accession" >> "$log_file_path"
        continue
    fi

    # Extract downloaded files
    echo "Extracting files..."
    if ! unzip -q "$download_zip_file" -d "$accession_temp_directory"; then
        echo "ERROR: Failed to extract $download_zip_file"
        echo "$(date '+%Y-%m-%d %H:%M:%S') - FAILED (extraction): $accession" >> "$log_file_path"
        continue
    fi

    # Remove zip file to save space
    rm "$download_zip_file"

    echo "Extraction completed"
    echo "Files extracted to: $accession_temp_directory"

    # Verify ncbi_dataset directory exists
    ncbi_data_directory="${accession_temp_directory}/ncbi_dataset/data"
    if [[ ! -d "$ncbi_data_directory" ]]; then
        echo "ERROR: Expected data directory not found: $ncbi_data_directory"
        echo "$(date '+%Y-%m-%d %H:%M:%S') - FAILED (no data dir): $accession" >> "$log_file_path"
        continue
    fi

    echo "Data directory verified: $ncbi_data_directory"
    echo "Download and extraction complete for $accession"

    echo ""
    echo "Extracting metadata..."

    # Find the assembly metadata file
    assembly_metadata_file=$(find "$ncbi_data_directory" -name "$expected_assembly_metadata" | head -n 1)

    if [[ ! -f "$assembly_metadata_file" ]]; then
        echo "ERROR: Assembly metadata file not found"
        echo "$(date '+%Y-%m-%d %H:%M:%S') - FAILED (no metadata): $accession" >> "$log_file_path"
        continue
    fi

    echo "Found metadata file: $assembly_metadata_file"

    # Extract organism name from JSON
    # Format: "organism":{"organismName":"Saccharomyces cerevisiae S288C",...}
    organism_name_raw=$(grep -o '"organismName":"[^"]*"' "$assembly_metadata_file" | head -n 1 | sed 's/"organismName":"//;s/"$//')

    if [[ -z "$organism_name_raw" ]]; then
        echo "ERROR: Could not extract organism name from metadata"
        echo "$(date '+%Y-%m-%d %H:%M:%S') - FAILED (no organism name): $accession" >> "$log_file_path"
        continue
    fi

    echo "Raw organism name: $organism_name_raw"

    # Remove spaces from organism name
    organism_name_clean=$(echo "$organism_name_raw" | sed 's/ /_/g')
    echo "Cleaned organism name: $organism_name_clean"

    # Extract strain from JSON
    # Format: "infraspecificNames":{"strain":"S288C"}
    strain_name=$(grep -o '"strain":"[^"]*"' "$assembly_metadata_file" | head -n 1 | sed 's/"strain":"//;s/"$//')

    # If no strain found, use NA
    if [[ -z "$strain_name" ]]; then
        strain_name="NA"
        echo "No strain found, using: $strain_name"
    else
        echo "Strain: $strain_name"
    fi

    # Extract assembly name for metadata file
    assembly_name=$(grep -o '"assemblyName":"[^"]*"' "$assembly_metadata_file" | head -n 1 | sed 's/"assemblyName":"//;s/"$//')
    echo "Assembly name: $assembly_name"

    # Extract tax ID for metadata file
    tax_id=$(grep -o '"taxId":[0-9]*' "$assembly_metadata_file" | head -n 1 | sed 's/"taxId"://')
    echo "Tax ID: $tax_id"

    # Sanitize accession - remove all dots and underscores
    accession_sanitized=$(echo "$accession" | sed 's/[._]//g')
    echo "Sanitized accession: $accession_sanitized"

    # Build filename prefix
    filename_prefix="${organism_name_clean}_${strain_name}_${accession_sanitized}"
    echo "Filename prefix: $filename_prefix"

    # Verify prefix looks reasonable
    if [[ ${#filename_prefix} -lt 10 ]]; then
        echo "ERROR: Filename prefix too short: $filename_prefix"
        echo "$(date '+%Y-%m-%d %H:%M:%S') - FAILED (bad prefix): $accession" >> "$log_file_path"
        continue
    fi

    echo "Metadata extraction complete"

    echo ""
    echo "Renaming and organizing files..."

    # Move all files from nested ncbi_dataset/data directory to temp root
    echo "Moving files from NCBI data directory..."
    find "$ncbi_data_directory" -type f -exec mv {} "$accession_temp_directory/" \;

    # Remove now-empty ncbi_dataset directory structure
    rm -rf "${accession_temp_directory}/ncbi_dataset"

    # Find and rename genome file
    echo "Renaming genome file..."
    genome_file_found=$(find "$accession_temp_directory" -name "*.fna" ! -name "*cds*" ! -name "*rna*" | head -n 1)

    if [[ -f "$genome_file_found" ]]; then
        genome_file_new="${accession_temp_directory}/${filename_prefix}_genome.fna"
        mv "$genome_file_found" "$genome_file_new"
        echo "  Genome: $(basename "$genome_file_new")"
    else
        echo "  WARNING: Genome file not found"
    fi

    # Find and rename CDS file
    echo "Renaming CDS file..."
    cds_file_found=$(find "$accession_temp_directory" -name "*cds*.fna" | head -n 1)

    if [[ -f "$cds_file_found" ]]; then
        cds_file_new="${accession_temp_directory}/${filename_prefix}_cds.fna"
        mv "$cds_file_found" "$cds_file_new"
        echo "  CDS: $(basename "$cds_file_new")"
    else
        echo "  WARNING: CDS file not found"
    fi

    # Find and rename protein file
    echo "Renaming protein file..."
    protein_file_found=$(find "$accession_temp_directory" -name "*.faa" | head -n 1)

    if [[ -f "$protein_file_found" ]]; then
        protein_file_new="${accession_temp_directory}/${filename_prefix}_protein.faa"
        mv "$protein_file_found" "$protein_file_new"
        echo "  Protein: $(basename "$protein_file_new")"
    else
        echo "  WARNING: Protein file not found"
    fi

    # Find and rename GFF annotation file
    echo "Renaming annotation file..."
    gff_file_found=$(find "$accession_temp_directory" -name "*.gff" -o -name "*.gff3" | head -n 1)

    if [[ -f "$gff_file_found" ]]; then
        # Keep original extension (.gff or .gff3)
        gff_extension="${gff_file_found##*.}"
        gff_file_new="${accession_temp_directory}/${filename_prefix}_annotation.${gff_extension}"
        mv "$gff_file_found" "$gff_file_new"
        echo "  Annotation: $(basename "$gff_file_new")"
    else
        echo "  WARNING: GFF annotation file not found"
    fi

    # Rename NCBI metadata files
    echo "Renaming metadata files..."

    if [[ -f "${accession_temp_directory}/${expected_assembly_metadata}" ]]; then
        mv "${accession_temp_directory}/${expected_assembly_metadata}" \
           "${accession_temp_directory}/${filename_prefix}_assembly_data_report.jsonl"
        echo "  Assembly metadata: ${filename_prefix}_assembly_data_report.jsonl"
    fi

    if [[ -f "${accession_temp_directory}/${expected_sequence_metadata}" ]]; then
        mv "${accession_temp_directory}/${expected_sequence_metadata}" \
           "${accession_temp_directory}/${filename_prefix}_sequence_report.jsonl"
        echo "  Sequence metadata: ${filename_prefix}_sequence_report.jsonl"
    fi

    if [[ -f "${accession_temp_directory}/${expected_dataset_catalog}" ]]; then
        mv "${accession_temp_directory}/${expected_dataset_catalog}" \
           "${accession_temp_directory}/${filename_prefix}_dataset_catalog.json"
        echo "  Dataset catalog: ${filename_prefix}_dataset_catalog.json"
    fi

    # Create simple human-readable metadata file
    echo "Creating metadata summary..."
    metadata_file_new="${accession_temp_directory}/${filename_prefix}_metadata.txt"

    cat > "$metadata_file_new" << EOF
        Organism: $organism_name_raw
        Strain: $strain_name
        Assembly: $assembly_name
        Accession: $accession
        Accession_Sanitized: $accession_sanitized
        Tax_ID: $tax_id
        Download_Date: $(date '+%Y-%m-%d %H:%M:%S')
        Download_Directory: $download_base_directory
        Filename_Prefix: $filename_prefix
EOF

    echo "  Metadata: $(basename "$metadata_file_new")"

    echo "File renaming complete"

    echo ""
    echo "Moving files to final location..."

    # Count files to move
    files_to_move_count=$(find "$accession_temp_directory" -type f | wc -l)
    echo "Files to move: $files_to_move_count"

    # Move all files from temp directory to final location
    find "$accession_temp_directory" -type f -exec mv {} "$download_base_directory/" \;

    echo "Files moved to: $download_base_directory"

    # Verify at least some files were moved
    moved_files_count=$(find "$download_base_directory" -type f -name "${filename_prefix}*" | wc -l)
    echo "Files with prefix '$filename_prefix' in final location: $moved_files_count"

    if [[ $moved_files_count -lt 3 ]]; then
        echo "WARNING: Expected more files to be moved (found only $moved_files_count)"
    fi

    # Clean up temporary directory for this accession
    echo "Cleaning up temporary directory..."
    rm -rf "$accession_temp_directory"
    echo "Removed: $accession_temp_directory"

    # Log successful completion
    echo "$(date '+%Y-%m-%d %H:%M:%S') - SUCCESS: $accession -> $filename_prefix" >> "$log_file_path"

    echo ""
    echo "Completed processing: $accession"
    echo "All files prefixed with: $filename_prefix"

done
# [After the for loop ends]

echo ""
echo "========================================"
echo "Download Summary"
echo "========================================"
echo "$(date '+%Y-%m-%d %H:%M:%S') - All processing complete" >> "$log_file_path"

# Count total files in final directory
total_files_count=$(find "$download_base_directory" -type f ! -path "*/temp_downloads/*" | wc -l)
echo "Total files in reference directory: $total_files_count"

# Show unique prefixes (one per genome)
echo ""
echo "Genomes processed:"
unique_prefixes=$(find "$download_base_directory" -type f -name "*_genome.fna" -exec basename {} \; | sed 's/_genome\.fna$//')

if [[ -z "$unique_prefixes" ]]; then
    echo "  WARNING: No genome files found"
else
    echo "$unique_prefixes" | while read prefix; do
        echo "  - $prefix"
    done
fi

# Clean up temp_downloads directory if empty
if [[ -d "$temp_download_directory" ]]; then
    remaining_temp_files=$(find "$temp_download_directory" -type f | wc -l)
    if [[ $remaining_temp_files -eq 0 ]]; then
        echo ""
        echo "Removing empty temporary download directory..."
        rmdir "$temp_download_directory"
        echo "Removed: $temp_download_directory"
    else
        echo ""
        echo "WARNING: Temporary directory still contains $remaining_temp_files files"
        echo "Location: $temp_download_directory"
    fi
fi


echo ""
echo "Log file location: $log_file_path"
echo "All operations completed successfully."
echo "========================================"
echo "All downloads completed"
echo "========================================"
