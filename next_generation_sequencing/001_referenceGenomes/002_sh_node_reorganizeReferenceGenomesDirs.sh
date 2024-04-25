#!/bin/bash
#DESCRIPTION: Reorganize the reference genomes to be in a single directory
#USAGE: Call script on directory with the reference genome directories obtained from ncbi
#ASSUMES: This script was written in conjuction with lab_utils/next_generation_sequencing/REFGENS
# Get the directory containing the downloaded files
while IFS= read -r dir; do
	download_dirs+=("$dir") # Replace with your actual directory
	echo $dir
done < <(find . -maxdepth 1 -type d | grep "/")
#("GCA_000146045.2" "GCA_000005845.2" "GCA_000001405.29" "GCA_002163515.1")

for download_dir in "${download_dirs[@]}"; do
	echo "Reorganizing $download_dir."
	
	# Extract organism name from assembly_data_report.jsonl (modified grep and sed)
	# SPECIFIC_SOLUTION: Some of the files extracted organismName more than once and I couldnt process the newline or space out. Workaround is to output the name to text file and grab first line.   
	organism_name=$(grep -m 1 -o '"organismName":"[^"]*' "$download_dir/ncbi_dataset/data/assembly_data_report.jsonl" | cut -d '"' -f4 | sed 's/ //g' )
	echo "$organism_name" > organism_name.txt
	organism_name=$(cat organism_name.txt | head -n 1)
	echo "Organism name is $organism_name" 
	rm organism_name.txt
	# Find all files recursively and move them to the root directory
	find "$download_dir/" -type f -exec mv -v {} "$download_dir/" \;
	
	#TEST:
	#find "$download_dir/" -type f -exec echo {} "$download_dir/" \;
	
	  # Rename fasta file based on organism name (optional)
	mv -v "$download_dir"/cds_from_genomic.fna "$download_dir"/cds.fna
	find "$download_dir" -name "*_genomic.fna" -exec mv -v {} "$download_dir"/"$organism_name"_refgenome.fna \;
	   
	echo "Files reorganized and renamed based on organism name."
	
	rm -rf "$download_dir"/ncbi_dataset
	mv "$download_dir" "$organism_name" && echo "Renamed dir and removed empty dir."
	
done

#awk '/^>/ {gsub(/chromosome/, "chr", $6); printf(">%s%s\n", $6, $7)}' data/REFGENS/SaccharomycescerevisiaeS288C/SaccharomycescerevisiaeS288C_refgenome.fna | sed 's/,.*//' > test.fasta
#awk '/^>/ {gsub(/chromosome/, "chr", $6); printf(">%s%s\n", $6, $7)} !/^>/ {print $0}' data/REFGENS/SaccharomycescerevisiaeS288C/SaccharomycescerevisiaeS288C_refgenome.fna | sed 's/,.*//' | head -n 10
