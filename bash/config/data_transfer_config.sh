#!/bin/bash

# Add to existing or create new configuration
declare -A BMC_CONFIG=(
    ["BASE_PATH"]="/net/%s/data/bmc/public/Bell/%s"
    ["LOCAL_BASE"]="$HOME/data"
    ["FASTQ_DIR"]="fastq"
)

declare -A RSYNC_CONFIG=(
    ["OPTIONS"]="-av"
    ["INCLUDES"]="--include '*/' --include '*.fastq'"
    ["EXCLUDES"]="--exclude '*'"
)

declare -A CLEANUP_PATTERNS=(
    ["DIRS"]=("*D24*" "infosite*")
    ["FILES"]=("*unmapped*.fastq")
)
