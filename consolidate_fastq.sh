#!/bin/bash

# Strict error handling
set -euo pipefail

# Function to validate FASTQ files
validate_fastq() {
    local file="$1"
    if ! [ -f "$file" ]; then
        echo "Error: $file not found"
        exit 1
    fi
    if ! [[ "$file" =~ \.fastq$ ]]; then
        echo "Error: $file is not a FASTQ file"
        exit 1
    fi
}

# Get unique IDs - now matching 5-6 digit sequences
unique_ids=$(ls *_sequence.fastq | grep -o '[0-9]\{5,6\}' | sort -u)

# Verify we found some IDs
if [ -z "$unique_ids" ]; then
    echo "Error: No valid IDs found in fastq files"
    exit 1
fi

echo "Found the following unique IDs:"
echo "$unique_ids"

# Process each unique ID
for id in $unique_ids; do
    echo "Processing ID: $id"
    
    # Find all files for this ID
    files=$(ls *"$id"*_sequence.fastq)
    
    # Validate all files before processing
    for file in $files; do
        validate_fastq "$file"
        echo "Validated: $file"
    done
    
    # Create output filename
    output_file="consolidated_${id}_sequence.fastq"
    
    # Consolidate files
    if cat $files > "$output_file"; then
        echo "Successfully created $output_file"
        
        # Verify the new file exists and has content
        if [ -s "$output_file" ]; then
            # Get size before and after
            original_size=$(du -b $files | awk '{sum += $1} END {print sum}')
            new_size=$(du -b "$output_file" | awk '{print $1}')
            
            echo "Original files total size: $original_size bytes"
            echo "New file size: $new_size bytes"
            
            if [ "$new_size" -gt 0 ]; then
                echo "Removing original files..."
                rm -f $files
                echo "Original files removed"
            else
                echo "Error: Consolidated file is empty"
                rm -f "$output_file"
                exit 1
            fi
        else
            echo "Error: Consolidated file is empty"
            rm -f "$output_file"
            exit 1
        fi
    else
        echo "Error during consolidation"
        rm -f "$output_file"
        exit 1
    fi
done

echo "All consolidation operations completed successfully"
