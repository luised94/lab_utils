#!/bin/bash
# bash/tests/core/test_error_conditions.sh
#
# Source test environment
source "$(dirname "${BASH_SOURCE[0]}")/test_setup.sh" || exit 1
setup_test_environment || exit 1

test_error_conditions() {
    echo "�� Testing error conditions"
    local failed=0
    
    # Test invalid log levels
    echo "�  �� Testing invalid log level"
    if log_message "INVALID" "Test" "/dev/null" 2>/dev/null; then
        echo "�  �  ? Should have failed"
        ((failed++))
    else
        echo "�  �  � Properly rejected"
    fi
    
    # Test protected paths
    echo "�  �� Testing protected paths"
    if acquire_lock "/etc/test" 2>/dev/null; then
        echo "�  �  ? Should have failed"
        ((failed++))
    else
        echo "�  �  � Properly rejected"
    fi
    
    # Test invalid configurations
    echo "�  �� Testing invalid config"
    local old_config="${CORE_CONFIG[LOG_LEVELS]}"
    CORE_CONFIG[LOG_LEVELS]=""
    if log_message "INFO" "Test" "/dev/null" 2>/dev/null; then
        echo "�  �  ? Should have failed"
        ((failed++))
    else
        echo "�  �  � Properly rejected"
    fi
    CORE_CONFIG[LOG_LEVELS]="$old_config"
    
    return $failed
}

[[ "${BASH_SOURCE[0]}" == "${0}" ]] && test_error_conditions
