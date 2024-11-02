#!/bin/bash
# bash/tests/core/run_tests.sh

#' Run Core Test Suite
run_core_tests() {

    local start_time=$(date +%s)
    echo "ษอออออออออออออออออออออออออออป"
    echo "บ Lab Utils Core Test Suite บ"
    echo "ศอออออออออออออออออออออออออออผ"
    
    local failed=0
    local test_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    
    # Basic verification first
    echo "รฤ Running basic verification"
    if ! "$test_dir/verify_initialize_lab_environment.sh"; then
        echo "ณ  ? Basic verification failed"
        echo "ภฤ Aborting further tests"
        return 1
    fi
    
    # Source test environment
    source "$test_dir/test_setup.sh" || {
        echo "? Failed to source test environment"
        return 1
    }
    setup_test_environment || {
        echo "? Failed to setup test environment"
        return 1
    }
    # Advanced tests
    local tests=(
        "test_logging_concurrent.sh"
        "test_lock_recovery.sh"
        "test_error_conditions.sh"
    )
    
    for test in "${tests[@]}"; do
        echo "รฤ Running $test"
        if ! "$test_dir/$test"; then
            ((failed++))
        fi
    done
    
    local duration=$(($(date +%s) - start_time))
    echo "Verification completed in ${duration}s"
    echo "ภฤ Test suite complete (Failed: $failed)"
    return $failed
}

[[ "${BASH_SOURCE[0]}" == "${0}" ]] && run_core_tests
