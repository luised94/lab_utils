#!/bin/bash
# Requires miniforge installation. See install_macs2_with_miniforge.sh
# Exit on error
set -e

# ANSI color codes
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# ASCII symbols for cross-environment compatibility
SUCCESS_MARK="[OK]"
FAILURE_MARK="[FAILED]"

# Function to check command status
check_command() {
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}${SUCCESS_MARK}${NC} $1"
    else
        echo -e "${RED}${FAILURE_MARK}${NC} $1"
        exit 1
    fi
}

echo "Starting MEME Suite 5.3.3 installation process..."

# Document system state
echo "System Information:" > ~/lab_utils/meme_install_log.txt
uname -a >> ~/lab_utils/meme_install_log.txt
gcc --version >> ~/lab_utils/meme_install_log.txt
python --version >> ~/lab_utils/meme_install_log.txt
perl --version >> ~/lab_utils/meme_install_log.txt

# Create conda environment
# 
conda create -n meme_env python=3.8 -y
check_command "Created conda environment"

# Activate environment
conda activate meme_env
check_command "Activated conda environment"

# Install system dependencies
echo "Installing system dependencies..."
sudo apt-get update
sudo apt-get install -y build-essential wget perl python zlib1g-dev libxml2-dev
check_command "Installed system dependencies"

# Document installed package versions
dpkg -l build-essential wget perl python zlib1g-dev libxml2-dev >> meme_install_log.txt

# Create installation directory
mkdir -p $HOME/meme_install
cd $HOME/meme_install
check_command "Created installation directory"

# Download MEME Suite
wget http://meme-suite.org/meme-software/5.3.3/meme-5.3.3.tar.gz
check_command "Downloaded MEME Suite"

# Extract archive
tar zxf meme-5.3.3.tar.gz
cd meme-5.3.3
check_command "Extracted MEME Suite"

# Configure installation
./configure --prefix=$HOME/meme \
    --with-url=http://meme-suite.org/ \
    --enable-build-libxml2 \
    --enable-build-libxslt
check_command "Configured MEME Suite"

# Build and install
make
check_command "Built MEME Suite"

make test
check_command "Tested MEME Suite"

make install
check_command "Installed MEME Suite"

# Set up environment
echo 'export PATH="$HOME/meme/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
check_command "Updated PATH"

# Verify installation
if meme -version >/dev/null 2>&1; then
    echo -e "${GREEN}${SUCCESS_MARK}${NC} MEME Suite installation verified"
else
    echo -e "${RED}${FAILURE_MARK}${NC} MEME Suite installation verification failed"
    exit 1
fi

# Export conda environment
conda env export > meme_env.yml
check_command "Exported conda environment"

# Clean up
cd $HOME
rm -rf $HOME/meme_install
check_command "Cleaned up installation files"

echo -e "\n${GREEN}Installation Complete${NC}"
echo "To activate MEME environment: conda activate meme_env"
echo "To deactivate: conda deactivate"
echo "Installation log saved in meme_install_log.txt"
