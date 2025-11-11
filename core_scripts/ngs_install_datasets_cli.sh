#!/bin/bash
#===============================================================================
# TITLE: Install NCBI Datasets Command Line Tools
# DESCRIPTION: Downloads and installs the NCBI datasets and dataformat CLI tools
# DATE: 2024-12-12
# VERSION: 2.0.0
#===============================================================================
set -euo pipefail


# Configuration
ncbi_base_url="https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64"
install_directory="/usr/local/bin"
tools_list=("datasets" "dataformat")

# Check if running with sudo
if [[ $EUID -ne 0 ]]; then
    echo "ERROR: This script must be run with sudo privileges"
    echo "Usage: sudo $0"
    exit 1
fi

# Check if datasets already installed
if command -v datasets &> /dev/null; then
    echo "NCBI datasets CLI is already installed"
    datasets --version
    echo ""
    read -p "Reinstall anyway? (y/n): " reinstall_confirm
    if [[ ! $reinstall_confirm =~ ^[Yy]$ ]]; then
        echo "Installation cancelled"
        exit 0
    fi
fi

# Check prerequisites
if ! command -v curl &> /dev/null; then
    echo "ERROR: curl is required but not installed"
    exit 1
fi

echo "========================================"
echo "Installing NCBI datasets CLI tools"
echo "========================================"
echo "Download URL: $ncbi_base_url"
echo "Install directory: $install_directory"
echo "Tools: ${tools_list[*]}"
echo ""

# Download and install each tool
for tool_name in "${tools_list[@]}"; do
    echo "Downloading $tool_name..."

    if ! curl -f -o "$tool_name" "$ncbi_base_url/$tool_name"; then
        echo "ERROR: Failed to download $tool_name"
        exit 1
    fi

    echo "Setting executable permissions..."
    chmod +x "$tool_name"

    echo "Moving to $install_directory..."
    mv "$tool_name" "$install_directory/"

    echo "Successfully installed: $tool_name"
    echo ""
done

# Verify installation
echo "========================================"
echo "Verifying Installation"
echo "========================================"

if datasets --version; then
    echo ""
    echo "Installation completed successfully!"
else
    echo "ERROR: Installation verification failed"
    exit 1
fi

# Reminder message
echo ""
echo "========================================"
echo "Usage Information"
echo "========================================"
cat << 'EOF'

The NCBI datasets CLI tools are now installed and ready to use.

To download reference genomes locally and transfer to cluster:
1. Run your genome download script on this machine
2. After download completes, rsync the data:

   rsync -avz --progress \
     ~/data/reference_genomes/ \
     YOUR_USER@cluster_hostname:~/data/reference_genomes/

For more information:
https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/

EOF
