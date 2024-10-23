#!/bin/bash

# File Management Configuration
declare -A FILE_TYPES=(
    ["SEQUENCE"]="fastq fq"
    ["ALIGNMENT"]="bam sam"
    ["VARIANT"]="vcf"
    ["INDEX"]="bai sai"
)

declare -A SEARCH_PATHS=(
    ["INCLUDE"]="*/code/* */script*/*"
    ["EXCLUDE"]="*/lib/* */renv/* */git/*"
)

declare -A DEFAULTS=(
    ["MAX_DEPTH"]=4
    ["BATCH_SIZE"]=25
    ["MIN_FREE_SPACE_GB"]=10
)
