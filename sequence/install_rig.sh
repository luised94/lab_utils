#!/bin/bash
# Purpose: Install rig R version manager.
# Usage: sudo ./install_rig.sh
# Date: 2025-11-25

# 1. Check if the script is running as root (sudo)
if [ "$EUID" -ne 0 ]; then
  echo "Error: This script must be run as root. Please use sudo."
  exit 1
fi

# 2. Check if 'rig' is already installed
if command -v rig &> /dev/null; then
  echo "Rig is already installed. Exiting."
  exit 0
fi

# 3. Check if 'curl' is installed
if ! command -v curl &> /dev/null; then
  echo "Error: curl is not installed. Please install curl first."
  exit 1
fi

# 4. Perform the installation
echo "Installing rig..."
curl -Ls "https://github.com/r-lib/rig/releases/download/latest/rig-linux-$(arch)-latest.tar.gz" \
  | tar xz -C /usr/local

# 5. Verify the installation by running rig
echo "Verifying installation..."
if rig --version &> /dev/null; then
  echo "Success! Rig is installed and running."
  # Optional: Print the version to the console
  rig --version
  echo "Use rig add to add a new R."
  echo "Example: rig add 4.4"
else
  echo "Error: Rig was extracted but failed to run."
  exit 1
fi
