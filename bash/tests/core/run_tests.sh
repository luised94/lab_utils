#!/bin/bash
# bash/tests/core/run_tests.sh

#' Run Core Test Suite
run_core_tests() {
    echo "浜様様様様様様様様様様様様様�"
    echo "� Lab Utils Core Test Suite �"
    echo "藩様様様様様様様様様様様様様�"
    
    local failed=0
    local test_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    
    # Basic verification first
    echo "団 Running basic verification"
    if ! "$test_dir/verify_initialize_lab_environment.sh"; then
        echo "�  ? Basic verification failed"
        echo "青 Aborting further tests"
        return 1
    fi
    
    # Advanced tests
    local tests=(
        "test_logging_concurrent.sh"
        "test_lock_recovery.sh"
        "test_error_conditions.sh"
    )
    
    for test in "${tests[@]}"; do
        echo "団 Running $test"
        if ! "$test_dir/$test"; then
            ((failed++))
        fi
    done
    
    echo "青 Test suite complete (Failed: $failed)"
    return $failed
}

[[ "${BASH_SOURCE[0]}" == "${0}" ]] && run_core_tests
