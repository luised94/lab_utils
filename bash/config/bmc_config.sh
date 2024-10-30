#!/bin/bash
declare -A BMC_CONFIG=(
    [SOURCE_FS]="/net/bmc-pub17/data/bmc/public/Bell"
    [TARGET_FS]="/net/bmc-pub14/data"
    [MIN_SPACE_GB]=50
)
## RSYNC settings
#[CLEANUP_DIRS]="*_fastqc *_stats"    # Add your cleanup patterns
#[CLEANUP_FILES]="*.html *.zip"        # Add your cleanup patterns
#[RSYNC_OPTIONS]="-av --progress"
#[RSYNC_INCLUDES]="--include '*/' --include '*.fastq'"
#[RSYNC_EXCLUDES]="--exclude '*'"
## Cleanup patterns
#[CLEANUP_DIRS]="*D24* infosite*"
#[CLEANUP_FILES]="*unmapped*.fastq"
## BMC-specific settings
##[BMC_BASE_PATH]="/net/%s/data/bmc/public/Bell/%s"
##[BMC_FASTQ_DIR]="fastq"
##[BMC_DEFAULT_SERVER]="bmc-pub17"
