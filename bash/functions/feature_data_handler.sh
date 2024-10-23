#!/bin/bash

source "../functions/git_handler.sh"

function setup_feature_directory() {
    local base_dir="${FEATURE_DATA[BASE_DIR]}"
    
    log_info "Setting up feature data directory: $base_dir"
    
    if [ ! -d "$base_dir" ]; then
        if ! mkdir -p "$base_dir"; then
            log_error "Failed to create directory: $base_dir"
            return 1
        }
        log_info "Created directory: $base_dir"
    }
    
    echo "$base_dir"
}

function download_rossi_data() {
    local base_dir="$1"
    local repo_config="${REPOSITORIES[ROSSI_2021]}"
    local target_dir="${base_dir}/${repo_config[dir]}"
    
    log_info "Downloading Rossi 2021 data"
    
    if ! clone_repository "${repo_config[url]}" "$target_dir" "${repo_config[depth]}"; then
        return 1
    }
    
    if ! verify_repository "$target_dir"; then
        return 1
    }
    
    log_info "Successfully downloaded Rossi 2021 data"
    return 0
}
