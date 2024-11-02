#!/bin/bash
# bash/tests/core/run_tests.sh

#' Run Core Test Suite
run_core_tests() {

    local start_time=$(date +%s)

    local repo_root
    if ! repo_root=$(git rev-parse --show-toplevel 2>/dev/null); then
        echo "? Not in a git repository"
        return 1
    fi
    
    local test_dir="$repo_root/bash/tests/core"
    local failed=0
    
    echo "� Lab Utils Core Test Suite �"
    ## Basic verification first
    #echo "�� Running basic verification"
    #if ! "$test_dir/verify_initialize_lab_environment.sh"; then
    #    echo "�  ? Basic verification failed"
    #    echo "�� Aborting further tests"
    #    return 1
    #fi
    
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
        echo "�� Running $test"
        if ! "$test_dir/$test"; then
            ((failed++))
            echo "�  ? Test failed"
        else
            echo "�  � Test passed"
        fi
    done
    
    local duration=$(($(date +%s) - start_time))
    echo "Verification completed in ${duration}s"
    echo "�� Test suite complete (Failed: $failed)"
    return $failed
}

[[ "${BASH_SOURCE[0]}" == "${0}" ]] && run_core_tests
