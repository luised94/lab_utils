#!/bin/bash

# Add to existing or create new configuration
declare -A FEATURE_DATA=(
    ["BASE_DIR"]="$HOME/data/feature_files"
    ["DEFAULT_DEPTH"]=1
)

declare -A REPOSITORIES=(
    ["ROSSI_2021"]={
        "url"="https://github.com/CEGRcode/2021-Rossi_Nature.git"
        "branch"="main"
        "depth"=1
        "dir"="rossi_2021"
    }
)
