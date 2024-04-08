#!/bin/bash

# Get the directory containing the downloaded files
download_dirs="./GCA_000146045.2_genome"  # Replace with your actual directory
#("GCA_000146045.2" "GCA_000005845.2" "GCA_000001405.29" "GCA_002163515.1")

for download_dir in "${download_dirs[@]}"; do
  echo "Reorganizing $download_dir."

  # Extract organism name from assembly_data_report.jsonl (modified grep and sed)
  organism_name=$(grep -o '"organismName":"[^"]*' "$download_dir/ncbi_dataset/data/assembly_data_report.jsonl" | cut -d '"' -f4 | sed 's/ //g')
  echo $organism_name

  # Find all files recursively and move them to the root directory
  find "$download_dir/" -type f -exec mv -v {} "$download_dir/" \;

  #Test: find "$download_dir/" -type f -exec echo {} "$download_dir/" \;

  # Rename fasta file based on organism name (optional)
  mv -v "$download_dir"/cds_from_genomic.fna "$download_dir"/cds.fna
  find "$download_dir" -name "*_genomic.fna" -exec mv -v {} "$download_dir"/"$organism_name"_refgenome.fna \;
   
  echo "Files reorganized and renamed based on organism name."

  rm -rf "$download_dir"/ncbi_dataset
  mv "$download_dir" "$organism_name" && echo "Renamed dir and removed empty dir."

done
