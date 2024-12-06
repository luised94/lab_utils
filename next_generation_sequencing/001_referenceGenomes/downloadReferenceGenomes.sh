#STATUS:  REMOVE.
#!/bin/bash

# Script Name: 001_downloadReferenceGenomes.sh
# Description: Download reference genomes to the REFGENS directory (relative path)
# Usage: ./001_downloadReferenceGenomes.sh

# Define base URL for NCBI genome downloads
base_url="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession"

# Define accessions to download (replace with your accessions)
accessions=("GCA_000146045.2" "GCA_000005845.2" "GCA_000001405.29" "GCA_002163515.1")

# Set a download directory (modify as needed)
download_dir="./REFGENS"

# Create the download directory if it doesn't exist
mkdir -p "$download_dir"

for accession in "${accessions[@]}"; do
    echo "Downloading genome for accession: $accession"

    # Construct download URL with desired annotation types
    download_url="${base_url}/${accession}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"

    # Download genome archive with error handling
    if curl -o "$download_dir/${accession}_genome.zip" "$download_url"; then
        echo "Download successful for $accession"
    else
        echo "Download failed for accession: $accession. Check network or accession validity."
        continue  # Skip to the next iteration if download fails
    fi

    # Unzip archive with suppressed output
    unzip > /dev/null 2>&1 "$download_dir/${accession}_genome.zip" -d "$download_dir/${accession}_genome" && rm "$download_dir/${accession}_genome.zip"

    echo "Download and extraction complete for accession: $accession"
done

echo "All requested genomes have been downloaded and extracted to $download_dir."
