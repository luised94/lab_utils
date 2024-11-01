#!/bin/bash

source "../config/environment_config.sh"

function validate_directory() {
    local dir="$1"
    
    log_info "Validating directory: $dir"
    
    if [ ! -d "$dir" ]; then
        log_error "Directory not found: $dir"
        return 1
    fi
    
    if [ ! -r "$dir" ]; then
        log_error "Directory not readable: $dir"
        return 1
    fi
    
    return 0
}

function find_function_files() {
    local base_dir="$1"
    local pattern="${FILE_PATTERNS[FUNCTIONS]}"
    
    log_info "Finding function files"
    
    local exclude_pattern=""
    for pattern in "${FILE_PATTERNS[EXCLUDE_PATTERNS][@]}"; do
        exclude_pattern+=" ! -name \"$pattern\""
    done
    
    eval "find \"$base_dir\" -type f -name \"$pattern\" $exclude_pattern"
}

function validate_function_file() {
    local file="$1"
    
    log_info "Validating function file: $file"
    
    if [ ! -f "$file" ]; then
        log_error "File not found: $file"
        return 1
    fi
    
    if [ ! -r "$file" ]; then
        log_error "File not readable: $file"
        return 1
    fi
    
    # Optional: Add syntax check
    if command -v bash > /dev/null; then
        if ! bash -n "$file"; then
            log_error "Syntax error in: $file"
            return 1
        fi
    fi
    
    return 0
}

function load_function_file() {
    local file="$1"
    
    log_info "Loading functions from: $file"
    
    if ! source "$file" 2>/dev/null; then
        log_error "Failed to source: $file"
        return 1
    fi
    
    return 0
}

function load_priority_functions() {
    local base_dir="$1"
    
    log_info "Loading priority functions"
    
    for file in "${LOAD_ORDER[PRIORITY][@]}"; do
        local full_path="$base_dir/$file"
        if [ -f "$full_path" ]; then
            load_function_file "$full_path" || continue
        else
            log_warning "Priority file not found: $file"
        fi
    done
}

function load_remaining_functions() {
    local base_dir="$1"
    
    log_info "Loading remaining functions"
    
    while IFS= read -r -d $'\0' file; do
        # Skip priority and optional files
        local basename=$(basename "$file")
        if [[ " ${LOAD_ORDER[PRIORITY][@]} ${LOAD_ORDER[OPTIONAL][@]} " =~ " $basename " ]]; then
            continue
        fi
        
        load_function_file "$file" || continue
        
    done < <(find_function_files "$base_dir")
}
