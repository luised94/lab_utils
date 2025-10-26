
#!/bin/bash
#===============================================================================
# TITLE: Install NCBI SRA-tools CLI
# DESCRIPTION: Downloads and configures SRA-tools for downloading SRA datasets
# AUTHOR: [Your Name]
# DATE: 2024-12-12
# VERSION: 1.0.0
#===============================================================================

#-------------------------------------------------------------------------------
# Configuration Variables
#-------------------------------------------------------------------------------
TOOLKIT_URL="https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz"
TOOLKIT_DIR=$(ls ${HOME} | grep sratoolkit | grep -v tar.gz)

#-------------------------------------------------------------------------------
# Check if toolkit is already installed
#-------------------------------------------------------------------------------
if command -v prefetch &> /dev/null && command -v fasterq-dump &> /dev/null; then
    echo "SRA-tools already installed and in PATH"
    exit 0
fi

if [[ -n "$TOOLKIT_DIR" ]]; then
    echo "SRA-tools found in ${HOME}/${TOOLKIT_DIR}"
    export PATH="${PATH}:${HOME}/${TOOLKIT_DIR}/bin"
    echo "Added SRA-tools to PATH for current session"
    echo "To make this permanent, add the following to your .bashrc:"
    echo "export PATH=\$PATH:\$HOME/${TOOLKIT_DIR}/bin"
    exit 0
fi

#-------------------------------------------------------------------------------
# Download and Install
#-------------------------------------------------------------------------------
echo "Downloading SRA toolkit..."
wget --output-document sratoolkit.tar.gz "$TOOLKIT_URL"

if [ $? -ne 0 ]; then
    echo "Error: Failed to download SRA toolkit"
    exit 1
fi

echo "Extracting SRA toolkit..."
tar -xzf sratoolkit.tar.gz -C "$HOME"

# Get the new toolkit directory name
TOOLKIT_DIR=$(ls ${HOME} | grep sratoolkit | grep -v tar.gz)

if [[ -z "$TOOLKIT_DIR" ]]; then
    echo "Error: Failed to extract SRA toolkit"
    exit 1
fi

# Add to PATH for current session
export PATH="${PATH}:${HOME}/${TOOLKIT_DIR}/bin"

echo "SRA toolkit installed successfully"
echo "Tools added to PATH for current session"
echo "To make this permanent, add the following to your .bashrc:"
echo "export PATH=\$PATH:\$HOME/${TOOLKIT_DIR}/bin"

# Cleanup
rm sratoolkit.tar.gz

# Test the installation
fastq-dump --stdout -X 2 SRR390728
