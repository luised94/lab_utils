#!/bin/bash
# bash/tests/core/test_logging_concurrent.sh

test_concurrent_logging() {
    local test_dir="/tmp/lab_utils_test_$$"
    local log_file="$test_dir/concurrent.log"
    local pids=()
    
    echo "ÃÄ Testing concurrent logging"
    mkdir -p "$test_dir"
    
    # Launch parallel processes
    for i in {1..5}; do
        (
            for j in {1..10}; do
                log_message "INFO" "Test $i-$j" "$log_file"
                sleep 0.1
            done
        ) &
        pids+=($!)
    done
    
    # Wait for completion
    for pid in "${pids[@]}"; do
        wait "$pid"
    done
    
    # Verify results
    local expected=50
    local actual=$(wc -l < "$log_file")
    
    if [[ "$actual" -eq "$expected" ]]; then
        echo "³  û All messages logged ($actual/$expected)"
        rm -rf "$test_dir"
        return 0
    else
        echo "³  ? Message count mismatch ($actual/$expected)"
        return 1
    fi
}

[[ "${BASH_SOURCE[0]}" == "${0}" ]] && test_concurrent_logging
