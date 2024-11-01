#!/bin/bash

declare -A SRA_CONFIG=(
    ["BASE_URL"]="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"
    ["DATA_DIR"]="$HOME/data"
)

# Study-specific configuration
declare -A EATON_2010=(
    ["BIOPROJECT"]="PRJNA117641"
    ["DESCRIPTION"]="ORC precisely positions nucleosomes at origins of replication"
    ["CONDITION"]="WT-G2-ORC"
    ["OUTPUT_FILE"]="nnNnH.fastq"
)

# Sample configuration
declare -A SAMPLES=(
    ["WT-G2-ORC-rep1.fastq.gz"]="SRR034475"
    ["WT-G2-ORC-rep2.fastq.gz"]="SRR034476"
)
