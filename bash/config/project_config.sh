#!/bin/bash
# bash/config/project_config.sh

declare -A PROJECT_CONFIG=(
    [REMOTE_HOST]="luria.mit.edu"
    [REMOTE_USER]="luised94"
    [REMOTE_PATH]="~/data"
    [REQUIRED_DIRS]="documentation fastq"
    [DEFAULT_LOG_ROOT]="$HOME/logs"
    [LOG_LEVELS]="TRACE DEBUG INFO WARNING ERROR FATAL"


    # BMC-specific settings
    [BMC_BASE_PATH]="/net/%s/data/bmc/public/Bell/%s"
    [BMC_FASTQ_DIR]="fastq"
    [BMC_DEFAULT_SERVER]="bmc-pub17"
    
    # RSYNC settings
    [RSYNC_OPTIONS]="-av --progress"
    [RSYNC_INCLUDES]="--include '*/' --include '*.fastq'"
    [RSYNC_EXCLUDES]="--exclude '*'"
    
    # Cleanup patterns
    [CLEANUP_DIRS]="*D24* infosite*"
    [CLEANUP_FILES]="*unmapped*.fastq"
    
    # File patterns
    [FASTQ_PATTERN]="*.fastq"
    [SAMPLE_PATTERN]="[0-9]{6}Bel"

    # Safety settings
    [REQUIRE_CONFIRMATION]="true"
    [STAGING_DIR_PREFIX]="bmc_staging"
)
