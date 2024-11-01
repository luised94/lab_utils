#!/bin/bash
# verify_core.sh

verify_core() {
    local start_time=$(date +%s)
    local failed=0
    
    # Test logging
    echo "Testing logging system..."
    local log_file
    log_file=$(initialize_logging "verify")
    
    # Test concurrent writes
    echo "Testing concurrent writes..."
    for i in {1..5}; do
        log_message "INFO" "Test message $i" "$log_file" &
    done
    wait
    
    # Test lock system
    echo "Testing lock system..."
    local lock_file="/tmp/test_$$.lock"
    if ! acquire_lock "$lock_file" 5; then
        echo "Lock acquisition failed"
        ((failed++))
    else
        release_lock "$lock_file"
    fi
    
    # Report results
    local duration=$(($(date +%s) - start_time))
    echo "Verification completed in ${duration}s"
    echo "Failed tests: $failed"
    
    return $failed
}

# Run verification
verify_core
