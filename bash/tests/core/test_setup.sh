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
    
    echo "�� Setting up test environment"
    echo "�  �� Repository root: $repo_root"
    
    # Source core configuration
    if ! source "$repo_root/bash/config/core_config.sh"; then
        echo "�  ? Failed to source core config"
        return 1
    fi
    
    # Source core modules
    local core_modules=(
        "logging.sh"
        "lock.sh"
    )
    
    for module in "${core_modules[@]}"; do
        if ! source "$repo_root/bash/core/$module"; then
            echo "�  ? Failed to source $module"
            return 1
        fi
    done
    
    # Verify environment
    echo "�  �� Verifying environment"
    echo "�  �  �� LOG_LEVELS: ${CORE_CONFIG[LOG_LEVELS]}"
    echo "�  �  �� MAX_MESSAGE_LENGTH: ${CORE_CONFIG[MAX_MESSAGE_LENGTH]}"
    
    return 0
}




























