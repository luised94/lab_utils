#!/bin/bash
# bash/tests/core/test_logging_concurrent.sh

test_concurrent_logging() {
    local test_dir="/tmp/lab_utils_test_$$"
    local log_file="$test_dir/concurrent.log"
    local pids=()
    
    # Initialize environment if running standalone
    if [[ -z "$LAB_UTILS_INITIALIZED" ]]; then
        source "$(dirname "${BASH_SOURCE[0]}")/test_setup.sh" || exit 1
        setup_test_environment || exit 1
    fi

    echo "[TEST] Testing concurrent logging"
    
    # Ensure test directory exists
    mkdir -p "$test_dir" || {
        echo "[FAIL] Failed to create test directory"
        return 1
    }
    
    # Touch log file to ensure it exists
    touch "$log_file" || {
        echo "[FAIL] Failed to create log file"
        return 1
    }

    # Export necessary functions for subprocesses
    export -f log_message
    export -f write_log_atomic
    
    # Launch parallel processes
    for i in {1..5}; do
        (
            # Source core config directly in subprocess
            source "${LAB_UTILS_ROOT}/bash/config/core_config.sh"
            
            for j in {1..10}; do
                log_message "INFO" "Test message $i-$j" "$log_file"
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
        echo "[PASS] All messages logged ($actual/$expected)"
        rm -rf "$test_dir"
        return 0
    else
        echo "[FAIL] Message count mismatch ($actual/$expected)"
        return 1
    fi
}

# Run test if executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    test_concurrent_logging
fi
