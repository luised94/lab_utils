#!/bin/bash

# Configuration
BASE_URL="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"
OUTPUT_DIR="${HOME}/data/100303Bel"
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
