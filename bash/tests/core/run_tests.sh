#!/bin/bash
# bash/tests/core/run_tests.sh

#' Run Core Test Suite
run_core_tests() {
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
    
    echo "ภฤ Test suite complete (Failed: $failed)"
    return $failed
}

[[ "${BASH_SOURCE[0]}" == "${0}" ]] && run_core_tests
