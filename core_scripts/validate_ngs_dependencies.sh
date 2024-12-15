#!/bin/bash

# Script name: validate_ngs_dependencies.sh
# Description: Validates dependencies for NGS analysis pipeline
# Created: 2024-12-15
# Author: Your Name

# ANSI color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# ASCII symbols for cross-environment compatibility
SUCCESS_MARK="[OK]"
FAILURE_MARK="[FAILED]"
WARNING_MARK="[WARNING]"

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to check if a module exists
module_exists() {
    module avail "$1" 2>&1 | grep -q "$1"
}

# Function to check R version
check_r_version() {
    local required_version="4.2.0"
    if command_exists R; then
        local installed_version=$(R --version | head -n1 | cut -d " " -f3)
        if [ "$(printf '%s\n' "$required_version" "$installed_version" | sort -V | head -n1)" = "$required_version" ]; then
            echo -e "${GREEN}${SUCCESS_MARK}${NC} R version $installed_version is compatible (required: $required_version)"
            return 0
        else
            echo -e "${RED}${FAILURE_MARK}${NC} R version $installed_version is not compatible (required: $required_version)"
            return 1
        fi
    else
        echo -e "${RED}${FAILURE_MARK}${NC} R is not installed"
        return 1
    fi
}

# Function to check Python version
check_python_version() {
    local required_version="2.7"
    if command_exists python; then
        local installed_version=$(python --version 2>&1 | cut -d " " -f2 | cut -d. -f1,2)
        if [ "$installed_version" = "$required_version" ]; then
            echo -e "${GREEN}${SUCCESS_MARK}${NC} Python version $installed_version is compatible"
            return 0
        else
            echo -e "${RED}${FAILURE_MARK}${NC} Python version $installed_version is not compatible (required: $required_version)"
            return 1
        fi
    else
        echo -e "${RED}${FAILURE_MARK}${NC} Python is not installed"
        return 1
    fi
}

# Function to check module version
check_module_version() {
    local module_name="$1"
    local required_version="$2"
    
    if module_exists "${module_name}/${required_version}"; then
        echo -e "${GREEN}${SUCCESS_MARK}${NC} ${module_name} version ${required_version} is available"
        return 0
    else
        echo -e "${RED}${FAILURE_MARK}${NC} ${module_name} version ${required_version} is not available"
        return 1
    fi
}

# Main validation process
echo "=== NGS Analysis Pipeline Dependency Validation ==="
echo "Checking required dependencies..."
echo

# Initialize error counter
errors=0

# Check R version
echo "Checking R..."
check_r_version || ((errors++))

# Check Python version
echo "Checking Python..."
check_python_version || ((errors++))

# Check NGS tools modules
echo "Checking NGS tools..."
check_module_version "bowtie2" "2.3.5.1" || ((errors++))
check_module_version "fastp" "0.20.0" || ((errors++))
check_module_version "fastqc" "0.11.5" || ((errors++))
check_module_version "deeptools" "3.0.1" || ((errors++))

# Special check for GATK (multiple versions might be acceptable)
echo "Checking GATK..."
if module_exists "gatk"; then
    echo -e "${GREEN}${SUCCESS_MARK}${NC} GATK is available"
else
    echo -e "${RED}${FAILURE_MARK}${NC} GATK is not available"
    ((errors++))
fi

# Check essential command line utilities
echo "Checking command line utilities..."
essential_utils=("awk" "sed" "grep" "cut" "sort" "uniq" "tr" "find" "xargs")
for util in "${essential_utils[@]}"; do
    if command_exists "$util"; then
        echo -e "${GREEN}${SUCCESS_MARK}${NC} $util is available"
    else
        echo -e "${RED}${FAILURE_MARK}${NC} $util is not available"
        ((errors++))
    fi
done

# Summary
echo
echo "=== Validation Summary ==="
if [ $errors -eq 0 ]; then
    echo -e "${GREEN}All dependencies are satisfied!${NC}"
    exit 0
else
    echo -e "${RED}Found $errors dependency issues that need to be resolved.${NC}"
    echo "Please install or load the missing dependencies before running the pipeline."
    exit 1
fi
