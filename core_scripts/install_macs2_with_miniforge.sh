#!/bin/bash
# Initial attempts using regular pip and python commands were unsuccesful.
# Switch to miniforge instead.
#
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

echo "Starting MACS2 installation process using Miniforge..."

# Download Miniforge
wget https://github.com/conda-forge/miniforge/releases/download/4.8.3-2/Miniforge3-4.8.3-2-Linux-x86_64.sh
check_command "Downloaded Miniforge installer"

# Install Miniforge with specified path and batch mode
bash Miniforge3-4.8.3-2-Linux-x86_64.sh -b -p ~/miniforge3
check_command "Installed Miniforge"

# Initialize conda in bash
~/miniforge3/bin/conda init bash
check_command "Initialized conda"

# Reload shell environment to recognize conda
source ~/.bashrc

# Verify conda installation
conda --version
check_command "Verified conda version"

# Create environment with Python 3.6
conda create -n macs2_env python=3.6 -y
check_command "Created conda environment"

# Activate the environment
conda activate macs2_env
check_command "Activated conda environment"

# Verify paths
echo "Verifying installation paths..."
which python
which pip
which conda

# Install MACS2
conda install -c bioconda macs2 -y
check_command "Installed MACS2"

# Verify MACS2 installation
macs2 --version
check_command "Verified MACS2 version"

# Test MACS2 functionality
if macs2 callpeak -h >/dev/null 2>&1; then
    echo -e "${GREEN}${SUCCESS_MARK}${NC} MACS2 help documentation accessible"
else
    echo -e "${RED}${FAILURE_MARK}${NC} MACS2 help documentation not accessible"
    exit 1
fi

# Clean up installer
rm Miniforge3-4.8.3-2-Linux-x86_64.sh
check_command "Cleaned up installation files"

echo -e "\n${GREEN}Installation Complete${NC}"
echo "To activate MACS2 environment: conda activate macs2_env"
echo "To deactivate: conda deactivate"
# After succesful installation and verification, you can capture the state with a yml file
# and reactivate easily. Although I have not tested this.
# conda env export > ~/lab_utils/macs2_env.yml
# conda env create -f ~/lab_utils/macs2_env.yml
# conda activate macs2_env
