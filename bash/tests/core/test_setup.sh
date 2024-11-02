# bash/tests/core/test_setup.sh

#!/bin/bash

#' Setup Test Environment
#' @description Source all required files for testing
setup_test_environment() {
    local test_root
    test_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
    
    # Source core files in correct order
    source "$test_root/core/initialize_lab_environment.sh" || {
        echo "? Failed to source initialization"
        return 1
    }
    
    echo "ÃÄ Test environment initialized"
    return 0
}
