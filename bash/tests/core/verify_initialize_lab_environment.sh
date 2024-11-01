#!/bin/bash
# verify_init.sh

# Quick initialization verification
verify_init() {
    local log_file
    
    # 1. Source initialization
    source "$HOME/lab_utils/bash/core/initialize_lab_environment.sh" || {
        echo "ERROR: Failed to source initialization script"
        return 1
    }
    
    # 2. Verify environment
    [[ -n "$LAB_UTILS_INITIALIZED" ]] || {
        echo "ERROR: Environment not initialized"
        return 1
    }
    
    # 3. Test logging
    log_file=$(initialize_logging "verify")
    [[ -f "$log_file" ]] || {
        echo "ERROR: Logging initialization failed"
        return 1
    }
    
    echo "û Initialization verified successfully"
    return 0
}

verify_init
