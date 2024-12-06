#!/bin/bash
declare -A BMC_CONFIG=(
    [SOURCE_FS]="/net/bmc-pub17/data/bmc/public/Bell"
    [TARGET_FS]="/net/bmc-pub14/data/bell/users/luised94"  # More specific
    [RSYNC_OPTIONS]="-av"
    [CLEANUP_DIRS]="*D24* infosite* *_fastqc *_stats"
    [CLEANUP_FILES]="*unmapped*.fastq *.html *.zip"
    [MIN_SPACE_GB]=50
    
)
## Cleanup patterns
## BMC-specific settings
##[BMC_BASE_PATH]="/net/%s/data/bmc/public/Bell/%s"
##[BMC_FASTQ_DIR]="fastq"
##[BMC_DEFAULT_SERVER]="bmc-pub17"
