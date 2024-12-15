#!/bin/bash

# Exit on error
set -e

# ANSI color codes
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# ASCII symbols for cross-environment compatibility
SUCCESS_MARK="[OK]"
FAILURE_MARK="[FAILED]"

echo "Starting MACS2 installation process..."

# Load required modules
if ! module load python3/3.6.4; then
    echo -e "${RED}${FAILURE_MARK}${NC} Failed to load Python module"
    exit 1
fi

if ! module load gnu/5.4.0; then
    echo -e "${RED}${FAILURE_MARK}${NC} Failed to load GNU module"
    exit 1
fi

# Create virtual environment
if python3 -m venv ~/macs2_env --clear --without-pip; then
    echo -e "${GREEN}${SUCCESS_MARK}${NC} Created virtual environment"
else
    echo -e "${RED}${FAILURE_MARK}${NC} Failed to create virtual environment"
    exit 1
fi

# Source the virtual environment
source ~/macs2_env/bin/activate
<<<<<<< HEAD
curl https://bootstrap.pypa.io/pip/3.6/get-pip.py -o get-pip.py



export PYTHONNOUSERSITE=1
export PIP_NO_CACHE_DIR=true
# Confirm virtual environment activation.
# (macs2_env) user@node dir$
which pip  # Should show ~/macs2_env/bin/pip
which python  # Should show ~/macs2_env/bin/python

# Upgrade pip and install numpy
echo "Installing dependencies..."
python3 -m pip install --upgrade pip
python -m pip install --ignore-installed numpy
pip install numpy

# Install MACS2
echo "Installing MACS2..."
pip install MACS2

# Verify installation
if macs2 --version; then
    echo -e "${GREEN}${SUCCESS_MARK}${NC} MACS2 installation completed successfully"
else
    echo -e "${RED}${FAILURE_MARK}${NC} MACS2 installation verification failed"
    exit 1
fi

echo "To activate the MACS2 environment, use: source ~/macs2_env/activate"

# Add alias to .bashrc (using more compatible syntax)
echo "Do you want to add an alias to activate MACS2 environment to your .bashrc? (y/n)"
read -r REPLY
if [ "$REPLY" = "y" ] || [ "$REPLY" = "Y" ]; then
    echo "alias activate_macs2='module load python3/3.6.4 gnu/5.4.0 && source ~/macs2_env/activate'" >> ~/.bashrc
    echo -e "${GREEN}${SUCCESS_MARK}${NC} Alias added. You can now use 'activate_macs2' to set up your MACS2 environment."
fi
