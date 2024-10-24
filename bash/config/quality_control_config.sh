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
