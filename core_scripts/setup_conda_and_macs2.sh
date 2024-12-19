#!/bin/bash

# conda_macs2_setup.sh
# Purpose: Initialize conda and set up MACS2 environment
# Usage: source conda_macs2_setup.sh

# Conda initialization block
__conda_setup="$('/home/luised94/miniforge3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/luised94/miniforge3/etc/profile.d/conda.sh" ]; then
        . "/home/luised94/miniforge3/etc/profile.d/conda.sh"
    else
        export PATH="/home/luised94/miniforge3/bin:$PATH"
    fi
fi
unset __conda_setup

# Function to verify MACS2 environment
setup_macs2_env() {
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

    echo "MACS2 environment successfully activated"
    return 0
}

# Run setup
setup_macs2_env
