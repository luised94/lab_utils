#!/bin/bash
# verify_init.sh

test_lock_management() {
    echo "ÃÄ Testing lock management"
    local test_lock="test_lock_$$"
    
    if acquire_lock "$test_lock" 5; then
        echo "³  û Lock acquired"
        if release_lock "$test_lock"; then
            echo "³  û Lock released"
            return 0
        fi
    fi
    echo "? Lock management failed"
    return 1
}

test_core_config() {
    local required_keys=(
        "VERSION"
        "LOG_LEVELS"
        "DEFAULT_LOG_ROOT"
        "LOCK_BASE_DIR"
    )
    
    for key in "${required_keys[@]}"; do
        [[ -n "${CORE_CONFIG[$key]}" ]] || return 1
    done
}

test_directory_permissions() {
    local dirs=(
        "${CORE_CONFIG[DEFAULT_LOG_ROOT]}"
        "${CORE_CONFIG[LOCK_BASE_DIR]}"
    )
    
    for dir in "${dirs[@]}"; do
        [[ -w "$dir" ]] || return 1
    done
}

test_guard_mechanism() {
    echo "ÃÄ Testing guard mechanism"
    (
        source "${LAB_UTILS_ROOT}/bash/core/initialize_lab_environment.sh"
        [[ -n "$LAB_UTILS_INITIALIZED" ]] || {
            echo "? Guard failed"
            return 1
        }
    )
    echo "³  û Guard verified"
    return 0
}
# Quick initialization verification
verify_init() {
    local failed=0
    echo "ÃÄ Starting verification"
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
    
    #Test directory permissions 
    test_directory_permissions || ((failed++))
    
    # 2. Verify environment
    [[ -n "$LAB_UTILS_INITIALIZED" ]] || {
        echo "ERROR: Environment not initialized"
        return 1
    }

    if ((failed > 0)); then
        echo "ÀÄ ? Verification failed ($failed errors)"
        return 1
    fi
    
    # Verify paths
    echo "ÃÄ Verifying paths"
    echo "³  LAB_UTILS_ROOT: $LAB_UTILS_ROOT"
    echo "³  Config directory: $LAB_UTILS_CONFIG_DIR"
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
    
    echo "û Initialization verified successfully"
    return 0
}

[[ "${BASH_SOURCE[0]}" == "${0}" ]] && verify_init
