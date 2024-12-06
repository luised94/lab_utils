# bash/tests/core/test_root_discovery.sh

#!/bin/bash

test_root_discovery() {
    echo "�� Testing root discovery"
    
    # Test from different directories
    local dirs=(
        "bash/core"
        "bash/tests/core"
        "R/core"
        "."
    )
    
    for dir in "${dirs[@]}"; do
        echo "�  �� Testing from: $dir"
        (
            cd "$dir" 2>/dev/null || {
                echo "�  �  ? Failed to change to directory"
                return 1
            }
            
            local root
            root="$(discover_lab_utils_root)" || {
                echo "�  �  ? Root discovery failed"
                return 1
            }
            
            echo "�  �  � Found root: $root"
        )
    done
    
    echo "�� Tests complete"
}

# Run if executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    source "../../core/initialize_lab_environment.sh"
    test_root_discovery
fi
