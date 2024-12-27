#!/bin/bash
# conda_macs2_setup.sh
# Purpose: Initialize conda and set up MACS2 environment
# Usage: source conda_macs2_setup.sh

# Conda initialization block with error handling
conda_init_failed=false
INSTALL_REFERENCE="$HOME/lab_utils/core_scripts/install_macs2.sh"

__conda_setup="$('/home/luised94/miniforge3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/luised94/miniforge3/etc/profile.d/conda.sh" ]; then
        . "/home/luised94/miniforge3/etc/profile.d/conda.sh"
    else
        if [ -f "/home/luised94/miniforge3/bin/conda" ]; then
            export PATH="/home/luised94/miniforge3/bin:$PATH"
        else
            echo "ERROR: Conda installation not found or incomplete"
            echo "Please refer to installation instructions in: $INSTALL_REFERENCE"
            conda_init_failed=true
        fi
    fi
fi
unset __conda_setup

if [ "$conda_init_failed" = true ]; then
    echo "ERROR: Failed to initialize conda environment"
    echo "You can set up MACS2 using: $INSTALL_REFERENCE"
    exit 1
fi

# Function to verify MACS2 environment
verify_macs2_env() {
    # Check if macs2_env exists
    if ! conda env list | grep -q "macs2_env"; then
        echo "ERROR: macs2_env not found" >&2
        return 1
    fi

    # Activate environment
    conda activate macs2_env

    # Verify MACS2 is available
    if ! command -v macs2 &> /dev/null; then
        echo "ERROR: MACS2 not found in environment" >&2
        return 1
    fi

    # Verify MACS2 is available and working (quietly)
    if ! macs2 --version &> /dev/null; then
        echo "ERROR: MACS2 not functioning properly" >&2
        return 1
    fi

    echo "MACS2 environment successfully activated"
    return 0
}

# Run setup
verify_macs2_env
