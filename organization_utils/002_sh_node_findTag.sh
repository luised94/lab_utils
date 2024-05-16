#!/bin/bash

# Enable strict mode 
set -euo pipefail

# Function to display usage
usage() {
	echo "Usage: $0 --tag <tag> [--directory <directory>] [--file-extensions <ext1,ext2,...>]"
	echo "Example: $0 --tag TODO --directory /path/to/search --file-extensions sh,R,py"
	exit 1
}

# Set default values
TAG="TODO"
DIRECTORY="."
FILE_EXTENSIONS=("sh" "R" "py")

# Parse command-line options
while [[ "$#" -gt 0 ]]; do
	case $1 in
		--tag)
			TAG="$2"
			shift 2
			;;
		--directory)
			DIRECTORY="$2"
			shift 2
			;;
		--file-extensions)
			IFS=',' read -r -a FILE_EXTENSIONS <<< "$2"
			shift 2
			;;
		-h|--help)
			usage
			;;
		*)
			echo "Unknown option: $1"
			usage
			;;
	esac
done
# Check if at least one argument is provided
if [ $# -lt 1 ]; then 
	echo "No tag provided. Defaulting to #TODO"
	TAG="TODO"
else 
	TAG=$1
fi 

# Set default directory to current if not provided
DIRECTORY="${2:-.}"

# Shift arguments to get file extensions
shift 2

# Default file extensions if none are provided
if [ $# -eq 0 ]; then 
	echo "No file extensions provided. Defaulting to sh, R and py."
	FILE_EXTENSIONS=("sh" "R" "py") 
else 
	FILE_EXTENSIONS=("$@")
fi 


# Build the find command with the provided or default file extensions
FIND_CMD="find \"${DIRECTORY}\" -type f "
for ext in "${FILE_EXTENSIONS[@]}"; do 
	FIND_CMD+="-o -name \"*.${ext}\" "
done

echo "Command to run $FIND_CMD"
# Use grep to search for #TODO at the start of a line within the found files 
# Handle spaces in filenames and ensure efficient processing

FIND_CMD=$(echo "$FIND_CMD" | sed 's/ -o / /')
echo "Final command to run:"
echo "eval "$FIND_CMD -print0" | xargs -0 grep -Hn "^#$TAG""

eval "$FIND_CMD -print0 | xargs -0 grep -Hn "^#$TAG""

#find "$DIRECTORY" -type f \( -name "*.sh" -o -name "*.R" \) -print0 | xargs -0 grep -Hn "^#$TAG"

#NUMBER_OF_INSTANCES=$(find "$DIRECTORY" -type f \( -name "*.sh" -o -name "*.R" \) -print0 | xargs -0 grep -Hn "^#$TAG" | wc -l)

#echo "$NUMBER_OF_INSTANCES instances of $TAG found"


