#!/bin/bash

declare -A NGS_CONFIG=(
    ["FILE_PREFIX"]="D24-"
    ["FILE_SUFFIX"]="_NA_sequence.fastq"
    ["EXCLUDE_PATTERNS"]="*unmapped* processed_*"
)

declare -A PATHS=(
    ["BASE_DATA"]="$HOME/data"
    ["FASTQ_SUBDIR"]="fastq"
    ["DOC_SUBDIR"]="documentation"
)

declare -A PATTERNS=(
    ["ID_PATTERN"]='[_-]'  # Pattern for ID extraction
    ["FASTQ_EXT"]=".fastq"
)
