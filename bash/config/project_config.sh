#!/bin/bash
# bash/config/project_config.sh

declare -A PROJECT_CONFIG=(
    [REMOTE_HOST]="luria.mit.edu"
    [REMOTE_USER]="luised94"
    [REMOTE_PATH]="~/data"
    [REQUIRED_DIRS]="documentation fastq"
    [DEFAULT_LOG_ROOT]="$HOME/logs"
    [LOG_LEVELS]="TRACE DEBUG INFO WARNING ERROR FATAL"




    # File patterns
    [FASTQ_PATTERN]="*.fastq"
    [SAMPLE_PATTERN]="[0-9]{6}Bel"

    # Safety settings
    [REQUIRE_CONFIRMATION]="true"
    [STAGING_DIR_PREFIX]="bmc_staging"
    #
    # FASTQ Processing
    [FASTP_BASE_PARAMS]="--cut_window_size 4 --cut_mean_quality 20 --n_base_limit 5 --average_qual 20 --qualified_quality_phred 20 --unqualified_percent_limit 50 --html /dev/null"
    [FASTP_STANDARD_PARAMS]="--length_required 50"
    [FASTP_EATON_PARAMS]="--length_required 20"
    
    # File Patterns
    [FASTQ_INPUT_PATTERN]="*.fastq"
    [FASTQ_EXCLUDE_PATTERNS]="*unmapped* processed_*"
    [FASTQ_OUTPUT_PREFIX]="processed_"

    # FASTQ Processing
    [FASTQ_PATTERN]="*.fastq"
    [BMC_FASTQ_ID_PATTERN]="[-_]"  # Pattern for splitting
    [BMC_ID_REGEX]="[0-9]{5,6}"    # Pattern for matching ID
    [FASTQ_EXCLUDE]="*unmapped* processed_*"
    [FASTQ_PREFIX]="processed_"
    [FASTQ_SUFFIX]=".fastq"

    # Directory Structure
    [FASTQ_DIR]="fastq"
    [PROCESSED_FASTQ_DIR]="processedFastq"
    [DOC_DIR]="documentation"
    
    # Module Requirements
    [REQUIRED_MODULES]="fastp/0.20.0"
)
