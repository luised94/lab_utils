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
    ["GENOME_PATTERNS"]=("GCF*.fna" "GCA*.fna")
)

# Add to existing genome configurations
declare -A GENOME_NAMING=(
    ["CDS_FILE"]="cds.fna"
    ["REFGENOME_SUFFIX"]="_refgenome.fna"
    ["ORIGINAL_CDS"]="cds_from_genomic.fna"
    ["GENOMIC_PATTERN"]="*_genomic.fna"
)

declare -A GENOME_PATHS=(
    ["NCBI_DATA"]="ncbi_dataset/data"
    ["ASSEMBLY_REPORT"]="assembly_data_report.jsonl"
)

declare -A CHROMOSOME_NAMING=(
    ["PATTERN_OLD"]="chromosome"
    ["PATTERN_NEW"]="chr"
    ["SEPARATOR"]=","
)
