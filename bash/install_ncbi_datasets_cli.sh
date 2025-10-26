#!/bin/bash
#===============================================================================
# TITLE: Install NCBI Datasets Command Line Tools
# DESCRIPTION: Downloads and installs the NCBI datasets and dataformat CLI tools
# AUTHOR: [Your Name]
# DATE: 2024-12-12
# VERSION: 1.0.0
#===============================================================================

#-------------------------------------------------------------------------------
# Configuration Variables
#-------------------------------------------------------------------------------
NCBI_BASE_URL="https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64"
INSTALL_DIR="/usr/local/bin"
TOOLS=(
    "datasets"
    "dataformat"
)

#-------------------------------------------------------------------------------
# Check Prerequisites
#-------------------------------------------------------------------------------
if ! command -v curl &> /dev/null; then
    echo "Error: curl is required but not installed. Please install curl first."
    exit 1
fi

if [[ $EUID -ne 0 ]] && [[ "$INSTALL_DIR" == "/usr/local/bin" ]]; then
    echo "Error: This script requires sudo privileges to install to $INSTALL_DIR"
    exit 1
fi

#-------------------------------------------------------------------------------
# Main Installation Process
#-------------------------------------------------------------------------------
echo "Starting NCBI datasets CLI tools installation..."

# Download and install each tool
for TOOL in "${TOOLS[@]}"; do
    echo "Downloading $TOOL..."
    if curl -f -o "$TOOL" "$NCBI_BASE_URL/$TOOL"; then
        echo "Successfully downloaded $TOOL"
    else
        echo "Error: Failed to download $TOOL"
        exit 1
    fi

    echo "Setting executable permissions for $TOOL..."
    chmod +x "$TOOL"

    echo "Moving $TOOL to $INSTALL_DIR..."
    if sudo mv "$TOOL" "$INSTALL_DIR/"; then
        echo "Successfully installed $TOOL to $INSTALL_DIR"
    else
        echo "Error: Failed to move $TOOL to $INSTALL_DIR"
        exit 1
    fi
done

#-------------------------------------------------------------------------------
# Verification
#-------------------------------------------------------------------------------
echo "Verifying installation..."
if datasets --version; then
    echo "datasets tool installed successfully"
else
    echo "Error: datasets tool installation verification failed"
    exit 1
fi

echo "Displaying datasets help information..."
datasets --help

echo "Installation completed successfully!"

#-------------------------------------------------------------------------------
# Usage Example
#-------------------------------------------------------------------------------
cat << 'EOF'
To use the NCBI datasets tools, try these commands:

datasets download genome help
datasets download taxonomy help
dataformat --help

For more information, visit:
https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/
EOF
