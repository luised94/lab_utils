#STATUS:
#!/usr/bin/env bash

set -o errexit
set -o nounset
set -o pipefail

TARGET_DIR="./files_to_remove/"  # User-specified target directory for moving files
SEARCH_DIR="."  # Define the directory to search in
MAX_DEPTH=4  # Define the maximum depth for the search

# Update file extensions for NGS data analysis
FILE_EXTENSIONS="*.fastq *.bam *.sam *.vcf *.fq"
SEARCH_PATHS="*/code/* */script*/*"

# Function to move files
move_files() {
    local target_dir=$1
    local search_dir=$2
    local max_depth=$3
    local -a conditions=()
    
    # Construct the find command conditions for NGS file extensions
    for ext in $FILE_EXTENSIONS; do
        echo "$ext"
        conditions+=("-o -name" "$ext")
    done
    
    # Remove the first '-o' from the conditions array
    conditions=("${conditions[@]:1}")
     
    for ext in $conditions; do
         echo "$ext"
    done


    # Execute the find command with constructed conditions
    find "$search_dir" -maxdepth "$max_depth" -type f \( "${conditions[@]}" \) \( -path "*/code/*" -o -path "*/script*/*" \) \
     | xargs -I {} echo {} | nl
     #   | head -n 25 \
     #  | xargs -I {} mv {} "$target_dir"
}

# Call the move_files function with the specified parameters
move_files "$TARGET_DIR" "$SEARCH_DIR" "$MAX_DEPTH"

find "." -type f \( -path "*/code/*" -o -path "*/script*/*" \) \( -name "*.fastq" -o -name "*.bam" -o -name "*.sam" -o -name "*.vcf" -o -name "*.fq" \)
find . -type d \( -path '*/lib/*' -o -path '*/renv/*' -o -path '*/git/*' \) -prune -o -type f -exec du -a {} + | sort -nr | head -n 1000 | awk -F/ '{print $NF}' | rev | cut -d. -f1 | rev | sort | uniq
