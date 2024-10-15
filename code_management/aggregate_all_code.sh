#!/bin/bash

set -euo pipefail

trap 'echo "Error on line $LINENO. Exit code: $?"; exit 1' ERR

REPO_FILE="repository_contents.xml"
VERBOSE=1

# Configurable directory exclusions
EXCLUDED_DIRS=(".git" "deprecatedCode" "renv")

# Function to echo only in verbose mode
verbose_echo() {
    if [[ $VERBOSE -eq 1 ]]; then
        echo "$@"
    fi
}

# Function to escape XML special characters
escape_xml() {
    sed -e 's/&/\&amp;/g' -e 's/</\&lt;/g' -e 's/>/\&gt;/g' -e "s/'/\&apos;/g" -e 's/"/\&quot;/g'
}

# Parse command line arguments
while getopts "q" opt; do
    case $opt in
        q) VERBOSE=0 ;;
        *) echo "Usage: $0 [-q]"; exit 1 ;;
    esac
done

verbose_echo "Starting repository code gathering and analysis with XML tagging..."

# Initialize XML file
echo '<?xml version="1.0" encoding="UTF-8"?>' > "$REPO_FILE"
echo '<repository>' >> "$REPO_FILE"

# Function to process files of a specific type
process_files() {
    local file_type=$1
    verbose_echo "Processing $file_type files..."
    
    # Check if git command is available
    if command -v git &> /dev/null && git rev-parse --is-inside-work-tree &> /dev/null; then
        root_directory=$(git rev-parse --show-toplevel)
    else
        verbose_echo "Git command not available or not in a Git repository."
        while true; do
            read -p "Do you want to assign the root_directory to '.'? (y/n) " answer
            case ${answer:0:1} in
                y|Y ) root_directory="."; break ;;
                n|N ) echo "Exiting."; exit 1 ;;
                * ) echo "Please answer yes or no." ;;
            esac
        done
    fi

    # Construct the find command with exclusions
    find_cmd="find \"$root_directory\" -type f -name \"*.$file_type\""
    for dir in "${EXCLUDED_DIRS[@]}"; do
        find_cmd+=" -not -path \"*/$dir/*\""
    done
    find_cmd+=" -print0"
    echo "$find_cmd"

    # Execute the find command
    while IFS= read -r -d '' file; do
        verbose_echo "Processing: $file"
        file_name=$(basename "$file")
        file_path=$(dirname "$file")
        mime_type=$(file -b --mime-type "$file")
        {
            echo "<file>"
            echo "  <name>$file_name</name>"
            echo "  <path>$file_path</path>"
            echo "  <type>$mime_type</type>"
            #echo "  <content>"
            echo "  <content><![CDATA["
            #cat "$file" || echo "Error: Unable to read $file"
            #echo "</content>"
            if [ -r "$file" ]; then
                cat "$file" | escape_xml
            else
                echo "Error: Unable to read $file" | escape_xml
            fi
            echo "]]></content>"
            echo "</file>"
        } >> "$REPO_FILE"
    done < <(eval "$find_cmd")
}

# Process files in the specified order, including Lua
for file_type in md sh R lua py; do
    process_files "$file_type"
done

# Close repository tag
echo '</repository>' >> "$REPO_FILE"

verbose_echo "Task completed. Check $REPO_FILE for gathered code."
