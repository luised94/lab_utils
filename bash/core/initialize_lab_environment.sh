#!/bin/bash
# bash/core/initialize_lab_environment.sh

#' Initialize Lab Utils Environment
#' @description Core initialization script for lab utilities
#' @export LAB_UTILS_ROOT, LAB_UTILS_INITIALIZED

# Guard against multiple inclusion
[[ -n "$LAB_UTILS_INITIALIZED" ]] && return

#' Discover Lab Utils Root Directory
#' @return String Absolute path to lab utils root
discover_lab_utils_root() {
    local script_path
    
    # Try readlink first (handles symlinks)
    if script_path="$(readlink -f "${BASH_SOURCE[0]}")"; then
        echo "$(cd "$(dirname "$script_path")/.." && pwd)"
    # Fallback to BASH_SOURCE
    elif [[ -n "${BASH_SOURCE[0]}" ]]; then
        echo "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
    # Last resort: use environment variable or home directory
    else
        echo "${LAB_UTILS_ROOT:-$HOME/lab_utils}"
    fi
}

# Set root directory
if [[ -z "$LAB_UTILS_ROOT" ]]; then
    LAB_UTILS_ROOT="$(discover_lab_utils_root)"
    export LAB_UTILS_ROOT
fi

# Verify critical directory exists
if [[ ! -d "$LAB_UTILS_ROOT" ]]; then
    echo "ERROR: Lab utils directory not found: $LAB_UTILS_ROOT" >&2
    return 1
fi

# Core configuration files (order matters)
readonly CORE_CONFIG_FILES=(
    "core_config.sh"      # Base configuration with logging and lock settings.
)

# Core module files (order matters)
readonly CORE_MODULES=(
    "logging.sh"          # Must be first
    "lock.sh"
    "path_utils.sh"
)

# Initialize core configuration
for config in "${CORE_CONFIG_FILES[@]}"; do
    config_path="${LAB_UTILS_ROOT}/bash/config/${config}"
    if [[ -f "$config_path" ]]; then
        source "$config_path"
    else
        echo "ERROR: Required configuration not found: ${config}" >&2
        echo "ERROR: Required configuration not found: ${LAB_UTILS_ROOT}" >&2
        return 1
    fi
done

# Initialize core modules
for module in "${CORE_MODULES[@]}"; do
    module_path="${LAB_UTILS_ROOT}/bash/core/${module}"
    if [[ -f "$module_path" ]]; then
        source "$module_path"
    else
        echo "ERROR: Required module not found: ${module}" >&2
        return 1
    fi
done

# Export common paths
export LAB_UTILS_CONFIG_DIR="${LAB_UTILS_ROOT}/bash/config"
export LAB_UTILS_CORE_DIR="${LAB_UTILS_ROOT}/bash/core"
export LAB_UTILS_MODULES_DIR="${LAB_UTILS_ROOT}/bash/modules"
export LAB_UTILS_SCRIPTS_DIR="${LAB_UTILS_ROOT}/bash/scripts"
export LAB_UTILS_TESTS_DIR="${LAB_UTILS_ROOT}/bash/tests"

# Initialize logging
if [[ -z "$LAB_UTILS_LOG_FILE" ]]; then
    initialize_logging "lab_utils_init"
fi

# Cleanup stale locks on initialization
#cleanup_stale_locks 2>/dev/null

# Mark as initialized
readonly LAB_UTILS_INITIALIZED=1

log_info "Lab environment initialized successfully"
