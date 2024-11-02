#!/bin/bash

# Test script for lock path validation
test_lock_validation() {
    local test_cases=(
        "my_process"                           # Should pass
        "../my_process"                        # Should fail
        "/tmp/lab_utils_locks/${USER}/test"   # Should pass
        "/etc/test"                           # Should fail
        "/tmp/other_user/test"                # Should fail
    )
    
    for test in "${test_cases[@]}"; do
        echo "Testing: $test"
        if validate_lock_path "$test"; then
            echo "û Passed"
        else
            echo "x Failed"
        fi
    done
}
