#!/bin/bash

# Enable strict mode 
set -euo pipefail

# Function to display usage
usage() {
	echo "Usage: $0 <tag> [directory]"
	echo "Example: $0 TODO /path/to/search"
	exit 1
}

# Check if at least one argument is provided
if [ $# -lt 1 ]; then 
	echo "No tag provided. Defaulting to #TODO"
	TAG="TODO"
else 
	TAG=$1
fi 

# Set default directory to current if not provided
DIRECTORY="${2:-.}"

# Set default directory to current if 
# Use find to locate .sh and .R files recursively
# Use grep to search for #TODO at the start of a line within the found files 
# Handle spaces in filenames and ensure efficient processing

find "$DIRECTORY" -type f \( -name "*.sh" -o -name "*.R" \) -print0 | xargs -0 grep -Hn "^#$TAG"

NUMBER_OF_INSTANCES=$(find "$DIRECTORY" -type f \( -name "*.sh" -o -name "*.R" \) -print0 | xargs -0 grep -Hn "^#$TAG" | wc -l)

echo "$NUMBER_OF_INSTANCES instances of $TAG found"


