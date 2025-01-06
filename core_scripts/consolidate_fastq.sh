#!/bin/bash
# Dependencies: Assumes fastq files where transfered to appropriate ~/data/<experiment_id> directory 

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

# Extract unique IDs using delimiter-based approach
# This specifically extracts the ID between the first and second dash after 'D24'
readarray -t unique_ids < <(ls *_sequence.fastq | awk -F'D24-' '{print $2}' | cut -d'-' -f1 | sort -u)

# Verify we found some IDs
if [ -z "$unique_ids" ]; then
    echo "Error: No valid IDs found in fastq files"
    exit 1
fi

echo "Found the following unique IDs:"
echo "----------------"
printf '%s\n' "${unique_ids[@]}" | xargs -n6 | sed 's/^/    /' | column -t
echo "----------------"
echo "Number of unique IDs: ${#unique_ids[@]}"
read -p "Proceed with job submission? (y/n): " confirm
confirm=$(echo "$confirm" | tr '[:upper:]' '[:lower:]')
if [[ "$confirm" != "y" ]]; then
    echo "Job submission cancelled"
    exit 0
fi

# Process each unique ID
for id in ${unique_ids[@]}; do
    echo "Processing ID: $id"
    # Find all files matching the specific pattern
    # Using more strict pattern matching to avoid false matches
    files=$(ls *"D24-${id}-"*_sequence.fastq 2>/dev/null || true)

    # Check if we found any files
    if [ -z "$files" ]; then
        echo "No files found for ID: $id"
        continue
    fi

    # Validate all files before processing
    for file in $files; do
        validate_fastq "$file"
        echo "Validated: $file"
    done

    # Create output filename
    output_file="consolidated_${id}_sequence.fastq"

    #echo "Successfully created $output_file"
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
