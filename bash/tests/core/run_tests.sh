#!/bin/bash
# bash/tests/core/run_tests.sh

run_core_tests() {
    local start_time=$(date +%s)
    local repo_root
    local failed=0

    if ! repo_root=$(git rev-parse --show-toplevel 2>/dev/null); then
        echo "[ERROR] Not in a git repository"
        return 1
    fi
    
    local test_dir="$repo_root/bash/tests/core"
    
    echo "[START] Lab Utils Core Test Suite"
    
    # Source and setup test environment
    #source "$test_dir/verify_initialize_lab_environment.sh" || {
    #    echo "[ERROR] Failed to source test environment"
    #    return 1
    #}
    
    # Source and setup test environment
    source "$test_dir/test_setup.sh" || {
        echo "[ERROR] Failed to source test environment"
        return 1
    }
    setup_test_environment || {
        echo "[ERROR] Failed to setup test environment"
        return 1
    }

    # Advanced tests
    local tests=(
        "test_logging_concurrent.sh"
        "test_lock_recovery.sh"
        "test_error_conditions.sh"
    )
    
    for test in "${tests[@]}"; do
        echo "[TEST] Running $test"
        if ! "$test_dir/$test"; then
            ((failed++))
            echo "[FAIL] | Test failed: $test"
        else
            echo "[PASS] | Test passed: $test"
        fi
    done
    
    local duration=$(($(date +%s) - start_time))
    echo "[INFO] Verification completed in ${duration}s"
    echo "[END] Test suite complete (Failed: $failed)"
    
    return $failed
}

# Run tests if executed directly
[[ "${BASH_SOURCE[0]}" == "${0}" ]] && run_core_tests
