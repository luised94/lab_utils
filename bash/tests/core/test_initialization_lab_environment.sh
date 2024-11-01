#!/bin/bash
# bash/tests/core/test_initialization.sh

# Source initialization script
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/bash/core/initialize_lab_environment.sh"

#' Test Lab Environment Initialization
#' @return Integer 0 if successful
test_initialization() {
    local test_log_dir="/tmp/lab_utils_test_$$"
    local failed=0
    
    echo "=== Testing Lab Environment Initialization ==="
    
    # Test 1: Guard Mechanism
    (
        source "${LAB_UTILS_ROOT}/bash/core/initialize_lab_environment.sh"
        [[ -n "$LAB_UTILS_INITIALIZED" ]] || exit 1
    ) || {
        echo "x Guard mechanism failed"
        ((failed++))
    }
    
    # Test 2: Path Resolution
    [[ -n "$LAB_UTILS_ROOT" ]] && [[ -d "$LAB_UTILS_ROOT" ]] || {
        echo "x Path resolution failed"
        ((failed++))
    }
    
    # Test 3: Configuration Loading
    [[ -n "${CORE_CONFIG[VERSION]}" ]] || {
        echo "x Configuration loading failed"
        ((failed++))
    }
    
    # Test 4: Logging Setup
    local test_log
    if test_log=$(initialize_logging "test_init"); then
        [[ -f "$test_log" ]] || {
            echo "x Log file creation failed"
            ((failed++))
        }
    else
        echo "x Logging initialization failed"
        ((failed++))
    fi
    
    # Test 5: Lock Management
    local test_lock="/tmp/test_lock_$$"
    if acquire_lock "$test_lock" 5; then
        release_lock "$test_lock" || {
            echo "x Lock release failed"
            ((failed++))
        }
    else
        echo "x Lock acquisition failed"
        ((failed++))
    fi
    
    # Cleanup
    rm -rf "$test_log_dir"
    
    # Report results
    if ((failed == 0)); then
        echo "û All initialization tests passed"
        return 0
    else
        echo "x $failed test(s) failed"
        return 1
    fi
}

# Run tests if executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    test_initialization
fi
