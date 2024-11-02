#!/bin/bash
source "$HOME/lab_utils/bash/config/bmc_config.sh"

remove_files_safely() {
    local pattern="$1"
    local log_file="$2"
    
    log_info "Removing items matching: $pattern" "$log_file"
    
    # Find matching items with sizes
    local items=$(find . -name "$pattern" -exec du -sh {} \; | \
                 sort -hr)
    
    if [[ -z "$items" ]]; then
        log_info "No items found matching: $pattern" "$log_file"
        return 0
    fi
    
    # Log items to be removed
    log_warning "Will remove:" "$log_file"
    echo "$items" | tee -a "$log_file"
    
    # Interactive confirmation if running interactively
    if [[ -t 0 ]]; then
        read -p "Proceed with deletion? (y/n) " -r
        [[ ! $REPLY =~ ^[Yy]$ ]] && return 0
    fi
    
    # Remove items
    echo "$items" | while read size item; do
        if rm -rf "${item#./}"; then
            log_info "Removed: $item ($size)" "$log_file"
        else
            log_error "Failed to remove: $item" "$log_file"
        fi
    done
}
