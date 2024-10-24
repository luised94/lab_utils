#!/bin/bash

# Add to existing FASTQ configurations
declare -A FASTP_PARAMS=(
    ["BASE_PARAMS"]=(
        "--cut_window_size 4"
        "--cut_mean_quality 20"
        "--n_base_limit 5"
        "--average_qual 20"
        "--qualified_quality_phred 20"
        "--unqualified_percent_limit 50"
        "--html /dev/null"
    )
    ["STANDARD"]=(
        "--length_required 50"
    )
    ["EATON"]=(
        "--length_required 20"
    )
)

declare -A FILE_PATTERNS=(
    ["INPUT"]="*.fastq"
    ["EXCLUDE_PATTERNS"]=("*unmapped*" "processed_*")
    ["OUTPUT_PREFIX"]="processed_"
)

declare -A OUTPUT_DIRS=(
    ["PROCESSED"]="processedFastq"
    ["LOGS"]="logs"
)

# Add to existing module configurations
declare -A REQUIRED_MODULES=(
    ["FASTP"]="fastp/0.20.0"
)
