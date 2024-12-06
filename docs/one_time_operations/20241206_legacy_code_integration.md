# Legacy Code Consolidation Process - December 2024

## Context
This document describes the process used to consolidate legacy code from 
`/mnt/c/Users/${WINDOWS_USER}/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff`
into our organized repository structure.

## Initial Setup
```bash
# Verify working directory
echo "Processing files from: /mnt/c/Users/${WINDOWS_USER}/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff"

# Create working branch
git checkout -b legacy_code_integration
mkdir -p ./docs/one_time_operations/20241206_legacy_code_integration.md

# Create necessary directories
mkdir -p code_review_holder/legacy/{R,Python,Shell,PyMOL}
mkdir -p data_sorting/to_organize

# Exclude system and package directories
find "/mnt/c/Users/${WINDOWS_USER}/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff" \
    -type f \
    ! -path "*/\.*" \
    ! -path "*/renv/*" \
    ! -path "*/lib/*" \
    ! -path "*/bin/*" \
    ! -path "*/node_modules/*" \
    ! -path "*/__pycache__/*" \
    | sed 's/.*\.//' | sort | uniq -c | sort -nr > file_extensions_inventory.txt
