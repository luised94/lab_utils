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
