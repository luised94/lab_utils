#!/bin/bash

# File patterns
SLURM_OUTPUT_PATTERN="slurm-*.out"
SLURM_ERROR_PATTERN="slurm-*.err"

# Directory settings
DEFAULT_SLURM_LOG_DIR="slurm_logs"
MAX_LOG_AGE_DAYS=30
MAX_DEPTH_SEARCH=2

# Size limits
MAX_FILE_SIZE_MB=100

# Batch settings
BATCH_SIZE=100

# Existing SLURM configuration, adding new parameters
declare -A SLURM_CONFIG=(
    ["NODES"]=1
    ["TASKS"]=1
    ["MEM_PER_CPU"]="50G"
    ["CPUS_PER_TASK"]=4
    ["EXCLUDE_NODES"]="c[5-22]"
    ["NICE"]=10000
    ["MAIL_TYPE"]="ALL"
    ["MAIL_USER"]="luised94@mit.edu"
)
#
# Add QC-specific configurations
declare -A QC_CONFIG=(
    ["OUTPUT_DIR"]="qualityControl"
    ["PROCESSED_PATTERN"]="processed_*.fastq"
    ["EXCLUDE_PATTERNS"]="*unmapped* processed_*"
)

# Extend existing module configuration
declare -A MODULES=(
    ["GNU"]="gnu/5.4.0"
    ["SAMTOOLS"]="samtools/1.10"
    ["FASTQC"]="fastqc/0.11.5"
)

# Add to existing SLURM configurations
declare -A SLURM_WRAPPER=(
    ["MAX_ARRAY_SIZE"]=16
    ["DEFAULT_ARRAY_FORMAT"]="1-%d%%16"
    ["SCRIPT_BASE_DIR"]="$HOME/lab_utils"
    ["DATA_BASE_DIR"]="$HOME/data"
)

declare -A LOG_PATTERNS=(
    ["OUTPUT"]="*.out"
    ["ERROR"]="*.err"
    ["SLURM"]="slurm-*.out"
)

declare -A TIME_FORMATS=(
    ["JOB_ID"]="%Y%m%d%M%S"
    ["LOG_SEARCH"]="%Y-%m-%d"
)
