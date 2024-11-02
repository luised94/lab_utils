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
