#!/bin/bash

set -euo pipefail

source "config/environment_config.sh"

function show_usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]
Initializes bash functions and environment

Options:
    -d, --dir DIR     Base directory (default: ${PATHS[BASE_DIR]})
    -v, --verbose     Enable verbose output
    -f, --force       Force reload all functions
    -h, --help        Show this help message
EOF
}

function initialize_environment() {
    local base_dir="$1"
    local functions_dir="$base_dir/${PATHS[FUNCTIONS_DIR]}"
    
    log_info "Initializing environment"
    log_info "Base directory: $base_dir"
    
    # Validate directories
    validate_directory "$functions_dir" || {
        log_error "Invalid functions directory"
        log_error "Please ensure lab_utils is properly installed"
        log_error "Use: git clone <repository_url>"
        return 1
    }
    
    # Load functions in order
    load_priority_functions "$functions_dir" || {
        log_error "Failed to load priority functions"
        return 1
    }
    
    load_remaining_functions "$functions_dir" || {
        log_warning "Some functions failed to load"
    }
    
    log_info "Environment initialization completed"
    return 0
}

function main() {
    local base_dir="${PATHS[BASE_DIR]}"
    local verbose=false
    local force=false
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -d|--dir) base_dir="$2"; shift 2 ;;
            -v|--verbose) verbose=true; shift ;;
            -f|--force) force=true; shift ;;
            -h|--help) show_usage; exit 0 ;;
            *) log_error "Unknown option: $1"; show_usage; exit 1 ;;
        esac
    done
    
    initialize_environment "$base_dir"
}

main "$@"
