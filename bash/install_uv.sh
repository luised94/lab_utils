#!/bin/bash
################################################################################
# Install uv
# Author: Luis | Date: 2025-11-12 | Version: 1.0.0
################################################################################
# PURPOSE:
#   Install uv as a python package and environment manager.
# USAGE:
#   chmod +x install_uv.sh;./install_uv.sh
# DEPENDENCIES:
#   Curl.
# OUTPUTS:
#   uv installed. Check with uv --version
################################################################################
if ! command -v curl &> /dev/null; then
    echo "Curl is required for this scripts."
    exit 0

fi

echo "Curl is installed..."

if command -v uv &> /dev/null; then
    echo "uv already installed ($(uv --version))"
    exit 0

fi

echo "uv is not installed. Proceeding with installation..."

curl -LsSf https://astral.sh/uv/install.sh | sh

if command -v uv &> /dev/null; then
    echo "uv installed successfully ($(uv --version))"
    exit 0
else
    echo "? Installation failed. Try: pip install uv"
    exit 1
fi
