#!/bin/bash

function organize_fastq_files() {
    local target_dir="$1"
    
    log_info "Organizing FASTQ files in: $target_dir"
    
    cd "$target_dir" || {
        log_error "Failed to change to directory: $target_dir"
        return 1
    }
    
    find . -type f -name "*.fastq" -exec mv {} . \; || {
        log_error "Failed to move FASTQ files"
        return 1
    }
    
    return 0
}

function safe_remove() {
    local type="$1"  # 'dir' or 'file'
    local pattern="$2"
    local staging_dir="tmp/staging_deleter_$$"
    
    log_info "Searching for ${type}s matching: $pattern"
    
    local items
    if [ "$type" = "dir" ]; then
        items=$(find . -type d -name "$pattern")
    else
        items=$(ls $pattern 2>/dev/null)
    fi
    
    if [ -z "$items" ]; then
        log_info "No ${type}s found matching: $pattern"
        return 0
    }
    
    log_warning "Found ${type}s to remove:"
    echo "$items"
    
    read -p "DELETE these ${type}s? (y/n) " -n 1 -r
    echo
    
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        log_info "Deletion aborted"
        return 0
    }
    
    mkdir -p "$staging_dir" || {
        log_error "Failed to create staging directory"
        return 1
    }
    
    echo "$items" | while read item; do
        mv "$item" "$staging_dir/" || log_warning "Failed to move: $item"
    done
    
    rm -rf "$staging_dir"
    log_info "Deletion complete"
    return 0
}
