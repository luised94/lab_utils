#!/bin/bash

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

# ASCII symbols for cross-environment compatibility
SUCCESS_MARK="[OK]"
FAILURE_MARK="[FAILED]"

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to check if a module exists
module_exists() {
    module avail "$1" 2>&1 | grep -q "$1"
}

echo "Validating MACS2 installation requirements..."

# Check if module command exists
if command_exists module; then
    echo -e "${GREEN}${SUCCESS_MARK}${NC} Module system is available"
else
    echo -e "${RED}${FAILURE_MARK}${NC} Module system is not available. This is required for loading Python and GNU compiler."
    exit 1
fi

# Check for Python 3.6.4 module
if module_exists python3/3.6.4; then
    echo -e "${GREEN}${SUCCESS_MARK}${NC} Python 3.6.4 module is available"
else
    echo -e "${RED}${FAILURE_MARK}${NC} Python 3.6.4 module is not available. This version is required for MACS2 installation."
    exit 1
fi

# Check for GNU 5.4.0 module
if module_exists gnu/5.4.0; then
    echo -e "${GREEN}${SUCCESS_MARK}${NC} GNU 5.4.0 module is available"
else
    echo -e "${RED}${FAILURE_MARK}${NC} GNU 5.4.0 module is not available. This version is required for compilation."
    exit 1
fi

# Check if virtual environment can be created
python3 -m venv test_env >/dev/null 2>&1
if [ $? -eq 0 ]; then
    echo -e "${GREEN}${SUCCESS_MARK}${NC} Virtual environment can be created"
    rm -rf test_env
else
    echo -e "${RED}${FAILURE_MARK}${NC} Unable to create virtual environment. Please check your permissions."
    exit 1
fi

# Check if pip is available
if command_exists pip; then
    echo -e "${GREEN}${SUCCESS_MARK}${NC} pip is available"
else
    echo -e "${RED}${FAILURE_MARK}${NC} pip is not available. It's required for installing MACS2."
    exit 1
fi

echo -e "\n${GREEN}All requirements are met. You can proceed with MACS2 installation.${NC}"
