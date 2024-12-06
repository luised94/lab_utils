#!/bin/bash
# bash/core/config_export.sh

export_core_config() {
    local serialized=""
    for key in "${!CORE_CONFIG[@]}"; do
        serialized+="${key}=${CORE_CONFIG[$key]};"
    done
    export LAB_UTILS_CONFIG_SERIALIZED="${serialized}"
    
    # Export individual keys
    for key in "${!CORE_CONFIG[@]}"; do
        export "LAB_UTILS_${key}=${CORE_CONFIG[$key]}"
    done
}

import_core_config() {
    declare -g -A CORE_CONFIG
    while IFS='=' read -r key value; do
        [[ -n "$key" ]] && CORE_CONFIG[$key]="$value"
    done < <(echo "$LAB_UTILS_CONFIG_SERIALIZED" | tr ';' '\n')
}
