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

    # Verify paths
    echo "Verifying paths:"
    echo "LAB_UTILS_ROOT: $LAB_UTILS_ROOT"
    echo "Config directory: $LAB_UTILS_CONFIG_DIR"
    
    # Check critical directories exist
    for dir in "bash/config" "bash/core" "bash/modules"; do
        if [[ ! -d "$LAB_UTILS_ROOT/$dir" ]]; then
            echo "ERROR: Required directory not found: $LAB_UTILS_ROOT/$dir"
            exit 1
        fi
    done

    # Check critical files exist
    for file in "core_config.sh" "logging_config.sh" "lock_config.sh"; do
        if [[ ! -f "$LAB_UTILS_ROOT/bash/config/$file" ]]; then
            echo "ERROR: Required configuration not found: $file"
            exit 1
        fi
    done
    #
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
