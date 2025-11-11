#!/bin/bash

# Configuration
BASE_URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"
OUTPUT_DIR="${HOME}/data/100303Bel/fastq"
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

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Download and process each sample
for sample in "${SAMPLES[@]}"; do
    IFS=':' read -r sample_name accession gsm <<< "$sample"
    output_file="${OUTPUT_DIR}/${sample_name}.fastq.gz"
    url="${BASE_URL}${accession:0:6}/${accession}/${accession}.fastq.gz"
    echo "Downloading $sample_name (Accession: $accession, GSM: $gsm)"

    curl -f -L --create-dirs -o "$output_file" "$url"

    if [ $? -eq 0 ]; then
        echo "Successfully downloaded $sample_name"

    else
        echo "Failed to download $sample_name" >&2
        continue

    fi

done

echo "All downloads completed."

# Group files by GSM
declare -A gsm_groups
for sample in "${SAMPLES[@]}"; do
    IFS=':' read -r name srr gsm <<< "$sample"

    if [[ ! ${gsm_groups[$gsm]} ]]; then
        gsm_groups[$gsm]=""

    fi

    gsm_groups[$gsm]+="${OUTPUT_DIR}/${name}.fastq.gz "

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
    output_file="${OUTPUT_DIR}/consolidated_${id}_sequence.fastq.gz"

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
