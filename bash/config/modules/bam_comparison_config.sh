#!/bin/bash

# Add to existing or create new configuration
declare -A BAMCOMPARE_PARAMS=(
    ["BIN_SIZE"]=10
    ["NORMALIZATION"]="CPM"
    ["SCALE_METHOD"]="readCount"
    ["GENOME_SIZE"]=12157105
    ["MIN_MAPPING_QUALITY"]=20
    ["OPERATION"]="ratio"
    ["IGNORE_REGIONS"]=("chrXII")
    ["PSEUDOCOUNT"]=1
)

declare -A OUTPUT_DIRS=(
    ["BIGWIG"]="bigwig"
    ["LOGS"]="logs"
)

declare -A R_CONFIG=(
    ["SCRIPT_PATH"]="$HOME/lab_utils/next_generation_sequencing/004_bamProcessing/determineInputForArrayId.R"
)

# Add to existing module configurations
declare -A REQUIRED_MODULES=(
    ["GNU"]="gnu/5.4.0"
    ["PYTHON"]="python/2.7.13"
    ["DEEPTOOLS"]="deeptools"
    ["R"]="r/4.2.0"
)
