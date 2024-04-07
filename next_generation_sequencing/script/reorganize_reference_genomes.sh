#!/bin/bash

# Get the directory containing the downloaded files
download_dir="./GCA_000146045.2_genome"  # Replace with your actual directory
echo $download_dir

# Extract organism name from assembly_data_report.jsonl (modified grep and sed)
organism_name=$(grep -o '"organismName":"[^"]*' "$download_dir/ncbi_dataset/data/assembly_data_report.jsonl" | cut -d '"' -f4 | sed 's/ //g')
echo $organism_name

# Find all files recursively and move them to the root directory
find "$download_dir/" -type f -exec mv -v {} "$download_dir/" \;

#Test: find "$download_dir/" -type f -exec echo {} "$download_dir/" \;
#Rename genomic file based on organism name (using xargs)
#find "$download_dir" -name "*_genomic.gff" -exec sh -c 'mv "$0" "{}"'"$organism_name"'.gff' \;

# Rename fasta file based on organism name (optional)
# You can uncomment and modify this block if needed
mv -v cds_from_genomic.fna cds.fna
find "$download_dir" -name "*_genomic.fna" -exec echo {} "$organism_name".fna \;
 
echo "Files reorganized and renamed based on organism name."






















