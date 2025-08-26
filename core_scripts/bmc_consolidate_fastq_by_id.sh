#!/bin/bash
# Dependencies: Assumes fastq files where transfered to appropriate ~/data/<experiment_id> directory 
# Requires bmc cleanup. Script will detect and exit if not.
# Usage: ./consolidate_fastq.sh [experiment_id]

# Strict error handling
set -euo pipefail

# ---- Argument Handling ----
if [[ ! $# -eq 1 ]]; then
    echo "Usage: $0 [experiment_id]"
    echo "Consolidates FASTQ files into single per-sample files"
    echo "  - Without arguments: uses current directory (must be ~/data/######Bel/fastq)"
    echo "  - With experiment_id: uses ~/data/<experiment_id>/fastq"
    exit 1
fi

original_dir="$(pwd)"
target_dir=""

if [[ $# -eq 1 ]]; then
    # Validate provided experiment_id format
    if [[ ! "$1" =~ ^[0-9]{6}Bel$ ]]; then
        echo "ERROR: Invalid experiment_id format - must be 6 digits followed by 'Bel'" >&2
        echo "Example: 230504Bel" >&2
        exit 2
    fi

    target_dir="${HOME}/data/$1/fastq"
    if [[ ! -d "$target_dir" ]]; then
        echo "ERROR: Directory not found: $target_dir" >&2
        echo "Check experiment_id or directory structure" >&2
        exit 2
    fi
else
    # Validate current directory format
    target_dir="$(pwd)"
    if [[ ! "$target_dir" =~ ^${HOME}/data/[0-9]{6}Bel/fastq$ ]]; then
        echo "ERROR: Current directory must be ~/data/######Bel/fastq" >&2
        echo "Path detected: $target_dir" >&2
        exit 2
    fi
fi

# ---- Main Execution ----
echo "Using FASTQ directory: ${target_dir}"
cd "$target_dir" || { echo "ERROR: Failed to enter directory"; exit 1; }

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

echo "Confirmed directory has been cleaned from other bmc files."

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
    echo "ERROR: No valid IDs found in fastq files" >&2
    echo "Ensure the unique ids logic is correct and no updates have occured to names." >&2
    exit 3
fi

echo "Processing ${#unique_ids[@]} sample IDs"
echo "Found the following unique IDs:"
echo "----------------"
printf '%s\n' "${unique_ids[@]}" | xargs -n6 | sed 's/^/    /' | column -t
echo "----------------"
read -rp "Proceed with job submission? (y/n): " confirm
confirm=$(echo "$confirm" | tr '[:upper:]' '[:lower:]')

if [[ "$confirm" != "y" ]]; then
    echo "Job submission cancelled"
    exit 4
fi

echo "Job confirmed. Proceed with consolidation."

# Process each unique ID
for id in "${unique_ids[@]}"; do
    echo "--------------------"
    echo "Processing ID: $id"
    # Find files using array and glob pattern
    files=( *"${id}"*_sequence.fastq )

    if [ ${#files[@]} -eq 0 ]; then
        echo "ERROR: No files found for ID: $id"
        continue
    fi

    if [ ${#files[@]} -eq 2 ]; then
        echo -e "[ERROR]: Found ${#files[@]} for $id\nExpected 2"
        continue
    fi

    # Validate all files before processing
    for file in "${files[@]}"; do
      if ! [ -f "$file" ]; then
          echo "Error: $file not found"
          exit 1
      fi
      if ! [[ "$file" =~ \.fastq$ ]]; then
          echo "Error: $file is not a FASTQ file"
          exit 1
      fi
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
                exit 5
            fi

            # Only remove if verification passed
            echo "Removing original files..."
            rm -f -- "${files[@]}"
            echo "Original files removed"
        else
            echo "Error: Consolidated file is empty"
            rm -f "$output_file"
            exit 5
        fi
    else
        echo "Error during consolidation"
        rm -f "$tmp_file" "$output_file"
        exit 5
    fi
    echo "--------------------"
done

echo "Consolidation completed successfully in ${target_dir}"
# Return to original directory
cd "$original_dir"
