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

# Find extensions.
find "/mnt/c/Users/${WINDOWS_USER}/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff" \
    -type f \
    ! -path "*/\.*" \
    ! -path "*/renv/*" \
    ! -path "*/lib/*" \
    ! -path "*/bin/*" \
    ! -path "*/node_modules/*" \
    ! -path "*/__pycache__/*" \
    | sed 's/.*\.//' | sort | uniq -c | sort -nr > file_extensions_inventory.txt

# Find large files.
find "/mnt/c/Users/${WINDOWS_USER}/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff" \
    -type f \
    ! -path "*/\.*" \
    ! -path "*/renv/*" \
    ! -path "*/lib/*" \
    ! -path "*/bin/*" \
    ! -path "*/node_modules/*" \
    ! -path "*/__pycache__/*" \
    -size +10M > large_files.txt

# Find expected scripts or documentation files.
find "/mnt/c/Users/Luised94/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff" -type f $$ \
    -name "*.R" -o \
    -name "*.r" -o \
    -name "*.qmd" -o \
    -name "*.py" -o \
    -name "*.pml" -o \
    -name "*.sh" -o \
    -name "*.ipynb" -o \
    -name "*.Rmd" \
) -print > code_files_to_move.txt

# Find potential data files
find "/mnt/c/Users/Luised94/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff" -type f $$ \
    -name "*.csv" -o \
    -name "*.xlsx" -o \
    -name "*.txt" -o \
    -name "*.tsv" -o \
    -name "*.rds" \
) > data_files.txt

# Find potentially missed files
comm -23 \
    <(find "/mnt/c/Users/Luised94/Dropbox (MIT)/Lab/Projects/automate-the-boring-stuff" -type f | sort) \
    <(cat data_files.txt code_files_to_move.txt | sort) > missed_files.txt

# Analyze missed files
cat missed_files.txt | sed 's/.*\.//' | sort | uniq -c | sort -nr > missed_files_analysis.txt

# Move files by type
while IFS= read -r file; do
    case "${file##*.}" in
        "R"|"Rmd") dest="code_review_holder/legacy/R/" ;;
        "py"|"ipynb") dest="code_review_holder/legacy/Python/" ;;
        "pml") dest="code_review_holder/legacy/PyMOL/" ;;
        "sh") dest="code_review_holder/legacy/Shell/" ;;
    esac
    cp -v "$file" "$dest"
done < code_files_to_move.txt

# Quick content check of unknown files
while IFS= read -r f; do
    echo "=== $f ===" >> unknown_files_analysis.txt
    file "$f" >> unknown_files_analysis.txt
    head -n 1 "$f" 2>/dev/null >> unknown_files_analysis.txt || echo "Binary file" >> unknown_files_analysis.txt
    echo >> unknown_files_analysis.txt
done < missed_files.txt
