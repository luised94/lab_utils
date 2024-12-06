#!/bin/bash
# Example script usage

# Initialize lab environment
if [[ -z "$LAB_UTILS_ROOT" ]]; then
    # Try common locations
    for dir in \
        "$HOME/lab_utils" \
        "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)" \
        "${LAB_UTILS_ROOT:-}"; do
        if [[ -f "$dir/bash/core/initialize_lab_environment.sh" ]]; then
            source "$dir/bash/core/initialize_lab_environment.sh"
            break
        fi
    done
    
    # Check initialization
    if [[ -z "$LAB_UTILS_INITIALIZED" ]]; then
        echo "ERROR: Failed to initialize lab environment" >&2
        exit 1
    fi
fi

# Script logic
main() {
    local log_file
    log_file=$(initialize_logging "$(basename "$0")")
    
    # Your code here
    log_info "Script running" "$log_file"
}

main "$@"
