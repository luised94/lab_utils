#!/bin/bash

# Add to existing QC configurations
declare -A QC_PATHS=(
    ["BASE_DIR"]="$HOME/data"
    ["QC_SUBDIR"]="qualityControl"
    ["MAX_DEPTH"]=1
)

declare -A FILE_PATTERNS=(
    ["FASTQC_ZIP"]="*.zip"
)

declare -A OPERATION_DEFAULTS=(
    ["CONFIRM_TIMEOUT"]=30
    ["UNZIP_BATCH_SIZE"]=10
    ["PRESERVE_ZIP"]=true
)

declare -A BAM_QC_OUTPUTS=(
    ["FLAGSTAT"]="_bamFlagstat.txt"
    ["QUICKCHECK"]="_bamQuickcheck.txt"
    ["STATS"]="_bamStats.txt"
)

declare -A SAMTOOLS_PARAMS=(
    ["FLAGSTAT"]="-O tsv"
    ["STATS"]=""
)

declare -A QC_DIRS=(
    ["OUTPUT"]="qualityControl"
    ["LOGS"]="logs"
)

# Add to existing module configurations
declare -A REQUIRED_MODULES=(
    ["GNU"]="gnu/5.4.0"
    ["SAMTOOLS"]="samtools/1.10"
    ["FASTQC"]="fastqc/0.11.5"
)
