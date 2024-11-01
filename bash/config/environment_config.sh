#!/bin/bash

declare -A PATHS=(
    ["BASE_DIR"]="$HOME/lab_utils"
    ["FUNCTIONS_DIR"]="bash/functions"
    ["CONFIG_DIR"]="bash/config"
    ["SCRIPTS_DIR"]="bash/scripts"
)

declare -A FILE_PATTERNS=(
    ["FUNCTIONS"]="*.sh"
    ["CONFIG"]="*.sh"
    ["EXCLUDE_PATTERNS"]=(".*" "_*" "test_*")
)

declare -A LOAD_ORDER=(
    ["PRIORITY"]=(
        "logging.sh"
        "lock_utils.sh"
        h"
    )
    ["OPTIONAL"]=(
        "experimental.sh"
        "deprecated.sh"
    )
)
