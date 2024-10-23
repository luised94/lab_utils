#!/bin/bash

source "../config/data_sources_config.sh"

# Add to existing git handler or create new
function validate_git() {
    if ! command -v git &>/dev/null; then
        log_error "Git is not installed"
        return 1
    }
    return 0
}

function clone_repository() {
    local repo_url="$1"
    local target_dir="$2"
    local depth="${3:-1}"
    
    if [ -d "$target_dir" ]; then
        log_warning "Directory already exists: $target_dir"
        return 1
    }
    
    log_info "Cloning repository: $repo_url"
    if ! git clone --depth="$depth" "$repo_url" "$target_dir"; then
        log_error "Failed to clone repository: $repo_url"
        return 1
    }
    
    log_info "Successfully cloned to: $target_dir"
    return 0
}

function verify_repository() {
    local repo_dir="$1"
    
    if [ ! -d "$repo_dir/.git" ]; then
        log_error "Not a git repository: $repo_dir"
        return 1
    }
    
    if ! git -C "$repo_dir" rev-parse HEAD &>/dev/null; then
        log_error "Invalid git repository: $repo_dir"
        return 1
    }
    
    return 0
}
