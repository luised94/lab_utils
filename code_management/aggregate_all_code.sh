#!/bin/bash

set -euo pipefail

trap 'echo "Error on line $LINENO. Exit code: $?"; exit 1' ERR

REPO_FILE="repository_contents.xml"
STATS_FILE="repository_stats.txt"

echo "Starting repository code gathering and analysis with XML tagging..."

# Initialize XML file
echo '<?xml version="1.0" encoding="UTF-8"?>' > "$REPO_FILE"
echo '<repository>' >> "$REPO_FILE"

# Function to process files of a specific type
process_files() {
    local file_type=$1
    echo "Processing $file_type files..."
    find . -type f -name "*.$file_type" -print0 | while IFS= read -r -d '' file; do
        echo "Processing: $file"
        file_name=$(basename "$file")
        file_path=$(dirname "$file")
        mime_type=$(file -b --mime-type "$file")
        
        {
            echo "<file>"
            echo "  <name>$file_name</name>"
            echo "  <path>$file_path</path>"
            echo "  <type>$mime_type</type>"
            echo "  <content><![CDATA["
            cat "$file" || echo "Error: Unable to read $file"
            echo "]]></content>"
            echo "</file>"
        } >> "$REPO_FILE"
    done
}

# Process files in the specified order, including Lua
for file_type in md sh R lua py; do
    process_files "$file_type"
done

# Close repository tag
echo '</repository>' >> "$REPO_FILE"

# Gather statistics
echo "Analyzing repository..."
{
    echo "File counts:"
    for file_type in md sh R lua py; do
        count=$(grep -c "<name>.*\.$file_type</name>" "$REPO_FILE")
        echo "$file_type: $count"
    done

    echo "Total lines of code:"
    grep -v '<.*>' "$REPO_FILE" | wc -l

    echo "Top 10 most frequent words:"
    grep -v '<.*>' "$REPO_FILE" | tr '[:upper:]' '[:lower:]' | tr -cs '[:alpha:]' '\n' | sort | uniq -c | sort -rn | head -n 10

    echo "Largest files:"
    grep -n '</file>' "$REPO_FILE" | while IFS=':' read -r line_num _; do
        file_name=$(sed -n "$(($line_num-4))p" "$REPO_FILE" | sed -e 's/.*<name>\(.*\)<\/name>.*/\1/')
        file_size=$(sed -n "$(($line_num-1))p" "$REPO_FILE" | wc -c)
        echo "$file_size $file_name"
    done | sort -rn | head -n 5

    echo "Total final file size:"
    stat -f%z "$REPO_FILE" || echo "Error: Unable to get file size"
} > "$STATS_FILE"

echo "Task completed. Check $REPO_FILE for gathered code and $STATS_FILE for statistics."

# Display stats
cat "$STATS_FILE"
