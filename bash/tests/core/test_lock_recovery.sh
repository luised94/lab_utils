#!/bin/bash
# bash/tests/core/test_lock_recovery.sh
#
# Source test environment
source "$(dirname "${BASH_SOURCE[0]}")/test_setup.sh" || exit 1
setup_test_environment || exit 1

test_lock_recovery() {
    echo "�� Testing lock recovery"
    local test_dir="/tmp/lab_utils_test_$$"
    local lock_file="$test_dir/stale.lock"
    source "$HOME/lab_utils/bash/config/core_config.sh"
    
    # Create stale lock
    mkdir -p "$lock_file"
    echo "99999" > "$lock_file/pid"
    
    # Test scenarios
    echo "�  �� Testing stale lock cleanup"
    if cleanup_stale_locks; then
        echo "�  �  � Cleanup successful"
    else
        echo "�  �  ? Cleanup failed"
        return 1
    fi
    
    echo "�  �� Testing forced release"
    if release_lock "stale" true; then
        echo "�  �  � Force release successful"
    else
        echo "�  �  ? Force release failed"
        return 1
    fi
    
    rm -rf "$test_dir"
    return 0
}

[[ "${BASH_SOURCE[0]}" == "${0}" ]] && test_lock_recovery
