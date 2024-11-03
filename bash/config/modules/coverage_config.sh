#!/bin/bash

# Add to existing or create new configuration
declare -A COVERAGE_PARAMS=(
    ["BIN_SIZE"]=10
    ["NORMALIZATION"]="CPM"
    ["MIN_MAPPING_QUALITY"]=20
    ["IGNORE_DUPLICATES"]=true
)
#--normalizeUsing {RPKM,CPM,BPM,RPGC}
declare -A FILE_PATTERNS=(
    ["BAM_SUFFIX"]="S288C.bam"
    ["BIGWIG_SUFFIX"]="_indivNorm.bw"
)

declare -A OUTPUT_DIRS=(
    ["BIGWIG"]="bigwig"
    ["LOGS"]="logs"
)

# Add to existing module configurations
declare -A REQUIRED_MODULES=(
    ["GNU"]="gnu/5.4.0"
    ["PYTHON"]="python/2.7.13"
    ["DEEPTOOLS"]="deeptools"
)
