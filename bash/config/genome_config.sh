#!/bin/bash

declare -A GENOME_CONFIG=(
    ["BASE_DIR"]="$HOME/data/REFGENS"
    ["LOG_DIR"]="$HOME/data/REFGENS/logs"
    ["GENOME_PATTERN"]="*_refgenome.fna"
    ["INDEX_SUFFIX"]="_index"
)

declare -A SLURM_CONFIG=(
    ["NODES"]=1
    ["TASKS"]=1
    ["MEM_PER_CPU"]="20G"
    ["EXCLUDE_NODES"]="c[5-22]"
    ["MAIL_TYPE"]="ALL"
    ["MAIL_USER"]="luised94@mit.edu"
)

declare -A MODULES=(
    ["GNU"]="gnu/5.4.0"
    ["BOWTIE2"]="bowtie2/2.3.5.1"
)

# Add NCBI-specific configurations
declare -A NCBI_CONFIG=(
    ["DOWNLOAD_INCLUDES"]="genome,rna,cds,protein,gff3,gtf"
    ["ACCESSIONS"]=(
        "GCF_000146045.2"  # S. cerevisiae S288C
        "GCF_000001405.40" # Human
        "GCF_000005845.2"  # E. coli
        "GCA_002163515.1"  # S cerevisaie W303
    )
    ["GENOME_PATTERNS"]="GCF*.fna GCA*.fna"
)
