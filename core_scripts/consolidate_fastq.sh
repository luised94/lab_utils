#!/bin/bash
# Dependencies: Assumes fastq files where transfered to appropriate ~/data/<experiment_id> directory 
# Run from appropriate directory.

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

# Validate cleanup state with single consolidated check
if [ "$(find . -mindepth 1 -type d | wc -l)" -gt 0 ] || \
   [ "$(find . -type f ! -name "*.fastq" | wc -l)" -gt 0 ] || \
   [ "$(find . -maxdepth 1 -type f -name "*.fastq" | wc -l)" -eq 0 ]; then
    echo "ERROR: Directory not properly cleaned up. Please ensure:
- No subdirectories exist
- Only FASTQ files remain
- FASTQ files are in current directory" >&2
    echo "Run ~/lab_utils/core_scripts/cleanup_bmc_directory.sh."
    exit 1
fi

# Extract unique IDs using delimiter-based approach
# This specifically extracts the ID between the first and second dash after 'D24'

readarray -t unique_ids < <(
    for f in *_sequence.fastq; do
        # Split filename components
        IFS='_-' read -ra parts <<< "${f##*/}"
        # Extract ID from standardized position (3rd component)
        printf "%s\n" "${parts[2]}"
    done | sort -u
)

# Verify we found some IDs
if [ ${#unique_ids[@]} -eq 0 ]; then
    echo "Error: No valid IDs found in fastq files"
    exit 1
fi

echo "Number of unique IDs: ${#unique_ids[@]}"
echo "Found the following unique IDs:"
echo "----------------"
printf '%s\n' "${unique_ids[@]}" | xargs -n6 | sed 's/^/    /' | column -t
echo "----------------"
read -p "Proceed with job submission? (y/n): " confirm
confirm=$(echo "$confirm" | tr '[:upper:]' '[:lower:]')

if [[ "$confirm" != "y" ]]; then
    echo "Job submission cancelled"
    exit 0
fi

# Process each unique ID
for id in "${unique_ids[@]}"; do
    echo "Processing ID: $id"
    # Find files using array and glob pattern
    files=( *"${id}"*_sequence.fastq )

    if [ ${#files[@]} -eq 0 ]; then
        echo "No files found for ID: $id"
        continue
    fi

    # Validate all files before processing
    for file in "${files[@]}"; do
        validate_fastq "$file"
        echo "Validated: $file"
    done

    output_file="consolidated_${id}_sequence.fastq"
    tmp_file="${output_file}.tmp"

    # Consolidate files with atomic write
    if cat -- "${files[@]}" > "$tmp_file"; then
        mv "$tmp_file" "$output_file"
        echo "Successfully created $output_file"

        # Verify file content
        if [ -s "$output_file" ]; then
            # Accurate size calculation
            original_size=$(wc -c "${files[@]}" | awk '/total/ {print $1}')
            new_size=$(wc -c < "$output_file")

            echo "Original files total size: $original_size bytes"
            echo "New file size: $new_size bytes"

            # Calculate original data checksum
            orig_checksum=$(cat "${files[@]}" | md5sum | cut -d' ' -f1)
            new_checksum=$(md5sum "$output_file" | cut -d' ' -f1)
            
            if [[ "$orig_checksum" != "$new_checksum" ]]; then
                echo "Error: Consolidated file checksum mismatch!" >&2
                echo "Expected: $orig_checksum" >&2
                echo "Actual:   $new_checksum" >&2
                exit 1
            fi

            # Only remove if verification passed
            echo "Removing original files..."
            rm -f -- "${files[@]}"
            echo "Original files removed"
        else
            echo "Error: Consolidated file is empty"
            rm -f "$output_file"
            exit 1
        fi
    else
        echo "Error during consolidation"
        rm -f "$tmp_file" "$output_file"
        exit 1
    fi
done

echo "All consolidation operations completed successfully"
