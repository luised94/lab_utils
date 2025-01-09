
#!/bin/bash
# verify_init.sh

test_lock_management() {
    echo "ÃÄ Testing lock management"
    local test_lock="test_lock_$$"
    
    if acquire_lock "$test_lock" 5; then
        echo "³  û Lock acquired"
        if release_lock "$test_lock"; then
            echo "³  û Lock released"
            return 0
        fi
    fi
    echo "? Lock management failed"
    return 1
}

test_core_config() {
    local required_keys=(
        "VERSION"
        "LOG_LEVELS"
        "DEFAULT_LOG_ROOT"
        "LOCK_BASE_DIR"
    )
    
    for key in "${required_keys[@]}"; do
        [[ -n "${CORE_CONFIG[$key]}" ]] || return 1
    done
}

test_directory_permissions() {
    local dirs=(
        "${CORE_CONFIG[DEFAULT_LOG_ROOT]}"
        "${CORE_CONFIG[LOCK_BASE_DIR]}"
    )
    
    for dir in "${dirs[@]}"; do
        [[ -w "$dir" ]] || return 1
    done
}

test_guard_mechanism() {
    echo "ÃÄ Testing guard mechanism"
    (
        source "${LAB_UTILS_ROOT}/bash/core/initialize_lab_environment.sh"
        [[ -n "$LAB_UTILS_INITIALIZED" ]] || {
            echo "? Guard failed"
            return 1
        }
    )
    echo "³  û Guard verified"
    return 0
}
# Quick initialization verification
verify_init() {
    local failed=0
    echo "ÃÄ Starting verification"
    local log_file
    
    # Source initialization
    source "$HOME/lab_utils/bash/core/initialize_lab_environment.sh" || {
        echo "? Failed to source initialization script"
        return 1
    }
    
    # Test core configuration
    test_core_config || ((failed++))
    
    # Test guard mechanism
    test_guard_mechanism || ((failed++))
    
    # Test lock management
    test_lock_management || ((failed++))
    
    #Test directory permissions 
    test_directory_permissions || ((failed++))
    
    # 2. Verify environment
    [[ -n "$LAB_UTILS_INITIALIZED" ]] || {
        echo "ERROR: Environment not initialized"
        return 1
    }

    if ((failed > 0)); then
        echo "ÀÄ ? Verification failed ($failed errors)"
        return 1
    fi
    
    # Verify paths
    echo "ÃÄ Verifying paths"
    echo "³  LAB_UTILS_ROOT: $LAB_UTILS_ROOT"
    echo "³  Config directory: $LAB_UTILS_CONFIG_DIR"
    #
    # Check critical directories exist
    for dir in "bash/config" "bash/core" "bash/modules"; do
        if [[ ! -d "$LAB_UTILS_ROOT/$dir" ]]; then
            echo "ERROR: Required directory not found: $LAB_UTILS_ROOT/$dir"
            exit 1
        fi
    done

    # Check critical files exist
    for file in "core_config.sh" ; do
        if [[ ! -f "$LAB_UTILS_ROOT/bash/config/$file" ]]; then
            echo "ERROR: Required configuration not found: $file"
            exit 1
        fi
    done
    #
    # 3. Test logging
    log_file=$(initialize_logging "verify_initialize_lab_environment")
    [[ -f "$log_file" ]] || {
        echo "ERROR: Logging initialization failed"
        return 1
    }
    
    echo "û Initialization verified successfully"
    return 0
}

[[ "${BASH_SOURCE[0]}" == "${0}" ]] && verify_init

# bash/tests/core/test_setup.sh

#!/bin/bash

#' Setup Test Environment
#' @description Source all required files for testing
setup_test_environment() {
    # Use git to find repository root
    local repo_root
    if ! repo_root=$(git rev-parse --show-toplevel 2>/dev/null); then
        echo "? Not in a git repository"
        return 1
    fi
    
    echo "ÃÄ Setting up test environment"
    echo "³  ÃÄ Repository root: $repo_root"
    
    # Source core configuration
    if ! source "$repo_root/bash/config/core_config.sh"; then
        echo "³  ? Failed to source core config"
        return 1
    fi
    
    # Source core modules
    local core_modules=(
        "logging.sh"
        "lock.sh"
    )
    
    for module in "${core_modules[@]}"; do
        if ! source "$repo_root/bash/core/$module"; then
            echo "³  ? Failed to source $module"
            return 1
        fi
    done
    
    # Verify environment
    echo "³  ÃÄ Verifying environment"
    echo "³  ³  ÃÄ LOG_LEVELS: ${CORE_CONFIG[LOG_LEVELS]}"
    echo "³  ³  ÀÄ MAX_MESSAGE_LENGTH: ${CORE_CONFIG[MAX_MESSAGE_LENGTH]}"
    
    return 0
}

# bash/tests/core/test_root_discovery.sh

#!/bin/bash

test_root_discovery() {
    echo "ÃÄ Testing root discovery"
    
    # Test from different directories
    local dirs=(
        "bash/core"
        "bash/tests/core"
        "R/core"
        "."
    )
    
    for dir in "${dirs[@]}"; do
        echo "³  ÃÄ Testing from: $dir"
        (
            cd "$dir" 2>/dev/null || {
                echo "³  ³  ? Failed to change to directory"
                return 1
            }
            
            local root
            root="$(discover_lab_utils_root)" || {
                echo "³  ³  ? Root discovery failed"
                return 1
            }
            
            echo "³  ³  û Found root: $root"
        )
    done
    
    echo "ÀÄ Tests complete"
}

# Run if executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    source "../../core/initialize_lab_environment.sh"
    test_root_discovery
fi

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

#!/bin/bash
# bash/tests/core/test_lock_recovery.sh
#
# Source test environment
source "$(dirname "${BASH_SOURCE[0]}")/test_setup.sh" || exit 1
setup_test_environment || exit 1

test_lock_recovery() {
    echo "ÃÄ Testing lock recovery"
    local test_dir="/tmp/lab_utils_test_$$"
    local lock_file="$test_dir/stale.lock"
    source "$HOME/lab_utils/bash/config/core_config.sh"
    
    # Create stale lock
    mkdir -p "$lock_file"
    echo "99999" > "$lock_file/pid"
    
    # Test scenarios
    echo "³  ÃÄ Testing stale lock cleanup"
    if cleanup_stale_locks; then
        echo "³  ³  û Cleanup successful"
    else
        echo "³  ³  ? Cleanup failed"
        return 1
    fi
    
    echo "³  ÃÄ Testing forced release"
    if release_lock "stale" true; then
        echo "³  ³  û Force release successful"
    else
        echo "³  ³  ? Force release failed"
        return 1
    fi
    
    rm -rf "$test_dir"
    return 0
}

[[ "${BASH_SOURCE[0]}" == "${0}" ]] && test_lock_recovery

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

#!/bin/bash
# bash/tests/core/test_error_conditions.sh
#
# Source test environment
source "$(dirname "${BASH_SOURCE[0]}")/test_setup.sh" || exit 1
setup_test_environment || exit 1

test_error_conditions() {
    echo "ÃÄ Testing error conditions"
    local failed=0
    
    # Test invalid log levels
    echo "³  ÃÄ Testing invalid log level"
    if log_message "INVALID" "Test" "/dev/null" 2>/dev/null; then
        echo "³  ³  ? Should have failed"
        ((failed++))
    else
        echo "³  ³  û Properly rejected"
    fi
    
    # Test protected paths
    echo "³  ÃÄ Testing protected paths"
    if acquire_lock "/etc/test" 2>/dev/null; then
        echo "³  ³  ? Should have failed"
        ((failed++))
    else
        echo "³  ³  û Properly rejected"
    fi
    
    # Test invalid configurations
    echo "³  ÃÄ Testing invalid config"
    local old_config="${CORE_CONFIG[LOG_LEVELS]}"
    CORE_CONFIG[LOG_LEVELS]=""
    if log_message "INFO" "Test" "/dev/null" 2>/dev/null; then
        echo "³  ³  ? Should have failed"
        ((failed++))
    else
        echo "³  ³  û Properly rejected"
    fi
    CORE_CONFIG[LOG_LEVELS]="$old_config"
    
    return $failed
}

[[ "${BASH_SOURCE[0]}" == "${0}" ]] && test_error_conditions
#!/bin/bash
# bash/tests/core/test_config_export.sh

test_config_export() {
    local test_dir="$repo_root/bash/tests/core"
    
    # First ensure environment is initialized
    source "$test_dir/test_setup.sh" || {
        echo "[ERROR] Failed to source test environment"
        return 1
    }
    setup_test_environment || {
        echo "[ERROR] Failed to setup test environment"
        return 1
    }

    echo "[TEST] Testing configuration export"
    
    # Test 1: Verify serialized config exists
    if [[ -z "$LAB_UTILS_CONFIG_SERIALIZED" ]]; then
        echo "[FAIL] Missing serialized configuration"
        # Debug output
        echo "[DEBUG] Environment variables:"
        env | grep LAB_UTILS
        return 1
    fi
    
    # Test 2: Verify config values are exported
    if [[ -z "$LAB_UTILS_LOG_LEVELS" ]]; then
        echo "[FAIL] Missing exported LOG_LEVELS"
        return 1
    fi
    
    # Test 3: Verify import function
    (
        # Clear existing config
        unset CORE_CONFIG
        declare -A CORE_CONFIG
        
        # Import config
        import_core_config
        
        # Verify critical values
        if [[ -z "${CORE_CONFIG[LOG_LEVELS]}" ]]; then
            echo "[FAIL] Failed to import LOG_LEVELS"
            return 1
        fi
        
        # Verify specific values match
        if [[ "${CORE_CONFIG[LOG_LEVELS]}" != "$LAB_UTILS_LOG_LEVELS" ]]; then
            echo "[FAIL] Imported config doesn't match exported values"
            echo "[DEBUG] Imported: ${CORE_CONFIG[LOG_LEVELS]}"
            echo "[DEBUG] Exported: $LAB_UTILS_LOG_LEVELS"
            return 1
        fi
    ) || return 1
    
    echo "[PASS] Configuration export verified"
    return 0
}

# Only run directly if not being sourced
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    # Get repository root
    repo_root=$(git rev-parse --show-toplevel 2>/dev/null) || {
        echo "[ERROR] Not in a git repository"
        exit 1
    }
    
    test_config_export
fi

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
        "test_config_export.sh"
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
