#!/bin/bash
# verify_init.sh

test_lock_management() {
    echo "�� Testing lock management"
    local test_lock="test_lock_$$"
    
    if acquire_lock "$test_lock" 5; then
        echo "�  � Lock acquired"
        if release_lock "$test_lock"; then
            echo "�  � Lock released"
            return 0
        fi
    fi
    echo "? Lock management failed"
    return 1
}

test_core_config() {
    echo "�� Testing core configuration"
    
    [[ -n "${CORE_CONFIG[VERSION]}" ]] || {
        echo "? Missing version"
        return 1
    }
    
    [[ -n "${CORE_CONFIG[LOG_LEVELS]}" ]] || {
        echo "? Missing log levels"
        return 1
    }
    
    echo "�  � Configuration verified"
    return 0
}

test_guard_mechanism() {
    echo "�� Testing guard mechanism"
    (
        source "${LAB_UTILS_ROOT}/bash/core/initialize_lab_environment.sh"
        [[ -n "$LAB_UTILS_INITIALIZED" ]] || {
            echo "? Guard failed"
            return 1
        }
    )
    echo "�  � Guard verified"
    return 0
}
# Quick initialization verification
verify_init() {
    local failed=0
    echo "�� Starting verification"
    local log_file
    
    # Source initialization
    source "$HOME/lab_utils/bash/core/initialize_lab_environment.sh" || {
        echo "? Failed to source initialization script"
        return 1
    }
    
    # Test core configuration
    test_core_config || ((failed++))
    
    # Test guard mechanism
    test_guard_mechanism || ((failed++))
    
    # Test lock management
    test_lock_management || ((failed++))
    
    # 2. Verify environment
    [[ -n "$LAB_UTILS_INITIALIZED" ]] || {
        echo "ERROR: Environment not initialized"
        return 1
    }

    if ((failed > 0)); then
        echo "�� ? Verification failed ($failed errors)"
        return 1
    fi
    
    # Verify paths
    echo "�� Verifying paths"
    echo "�  LAB_UTILS_ROOT: $LAB_UTILS_ROOT"
    echo "�  Config directory: $LAB_UTILS_CONFIG_DIR"
    #
    # Check critical directories exist
    for dir in "bash/config" "bash/core" "bash/modules"; do
        if [[ ! -d "$LAB_UTILS_ROOT/$dir" ]]; then
            echo "ERROR: Required directory not found: $LAB_UTILS_ROOT/$dir"
            exit 1
        fi
    done

    # Check critical files exist
    for file in "core_config.sh" ; do
        if [[ ! -f "$LAB_UTILS_ROOT/bash/config/$file" ]]; then
            echo "ERROR: Required configuration not found: $file"
            exit 1
        fi
    done
    #
    # 3. Test logging
    log_file=$(initialize_logging "verify_initialize_lab_environment")
    [[ -f "$log_file" ]] || {
        echo "ERROR: Logging initialization failed"
        return 1
    }
    
    echo "� Initialization verified successfully"
    return 0
}

verify_init
