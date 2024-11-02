#!/bin/bash
# bash/tests/core/test_config_export.sh

test_config_export() {
    local test_dir="$repo_root/bash/tests/core"
    
    # First ensure environment is initialized
    source "$test_dir/test_setup.sh" || {
        echo "[ERROR] Failed to source test environment"
        return 1
    }
    setup_test_environment || {
        echo "[ERROR] Failed to setup test environment"
        return 1
    }

    echo "[TEST] Testing configuration export"
    
    # Test 1: Verify serialized config exists
    if [[ -z "$LAB_UTILS_CONFIG_SERIALIZED" ]]; then
        echo "[FAIL] Missing serialized configuration"
        # Debug output
        echo "[DEBUG] Environment variables:"
        env | grep LAB_UTILS
        return 1
    fi
    
    # Test 2: Verify config values are exported
    if [[ -z "$LAB_UTILS_LOG_LEVELS" ]]; then
        echo "[FAIL] Missing exported LOG_LEVELS"
        return 1
    fi
    
    # Test 3: Verify import function
    (
        # Clear existing config
        unset CORE_CONFIG
        declare -A CORE_CONFIG
        
        # Import config
        import_core_config
        
        # Verify critical values
        if [[ -z "${CORE_CONFIG[LOG_LEVELS]}" ]]; then
            echo "[FAIL] Failed to import LOG_LEVELS"
            return 1
        fi
        
        # Verify specific values match
        if [[ "${CORE_CONFIG[LOG_LEVELS]}" != "$LAB_UTILS_LOG_LEVELS" ]]; then
            echo "[FAIL] Imported config doesn't match exported values"
            echo "[DEBUG] Imported: ${CORE_CONFIG[LOG_LEVELS]}"
            echo "[DEBUG] Exported: $LAB_UTILS_LOG_LEVELS"
            return 1
        fi
    ) || return 1
    
    echo "[PASS] Configuration export verified"
    return 0
}

# Only run directly if not being sourced
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    # Get repository root
    repo_root=$(git rev-parse --show-toplevel 2>/dev/null) || {
        echo "[ERROR] Not in a git repository"
        exit 1
    }
    
    test_config_export
fi
