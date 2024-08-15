#!/bin/bash
# Description: This script searches for tags in the format #TAG on a given directory and given file extensions.
# Run 002_sh_node_findTag.sh -h for usage information

# Enable strict mode 
set -euo pipefail

#Find the readme file. Assumes it is in home directory
#OPTIMIZE: Not very robust but unclear how unless lab_utils is a unique directory that can be found.
README_WITH_TAGS=$(find "$HOME/lab_utils" -maxdepth 1 -type f -name "README.md")
#echo $README_WITH_TAGS

#Extract tags from README file, like mining gems from a document!
TAGS_IN_README=$(sed -n '/^## TAGS/,/^#/p' "$README_WITH_TAGS" | grep -v "#" | grep ":" | sed '/^[[:space:]]*$/d' | sed 's/:.*//' | tr '\n' ' ')
#echo $TAGS_IN_README

# Function to display usage
usage() {
	echo -e "Usage: $0 --tag <tag> [--directory <directory>] [--file-extensions <ext1,ext2,...>]\n"
	echo -e "Example: $0 --tag TODO --directory /path/to/search --file-extensions sh,R,py\n"
	echo "File extensions must be comma separated."
	echo -e "Options dont have to be in a particular order.\n"
	echo -e "Proper tags are: $TAGS_IN_README"
	exit 1
}


# Set default values
TAG="TODO"
DIRECTORY="."
FILE_EXTENSIONS=("sh" "R" "py")

# Parse command-line options
# As long as there are arguments, process using case argument. 
#NOTE shift helps process the arguments like a conveyor belt. 
while [[ "$#" -gt 0 ]]; do
	case $1 in
#TODO Use IFS and add for loop to search for multiple tags. 
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
			echo -e "Unknown option: $1 \n"
			usage
			;;
	esac
done

#If the tag is not in the TAGS_IN_README array, print message. Does not stop file for now
if [[ ! " ${TAGS_IN_README[*]} " =~ [[:space:]]${TAG}[[:space:]] ]]; then
	echo "$TAG is not in README tags, may not find result"
	echo "Proper tags are: $TAGS_IN_README"
fi

# Build the find command with the provided or default file extensions
FIND_CMD="find \"${DIRECTORY}\" -type f \( "
for ext in "${FILE_EXTENSIONS[@]}"; do 
	FIND_CMD+="-o -name \"*.${ext}\" "
done

#echo "Command to run $FIND_CMD"
# Use grep to search for #TODO at the start of a line within the found files 
# Handle spaces in filenames and ensure efficient processing

# Remove the initial -o from the find command
FIND_CMD=$(echo "$FIND_CMD \)" | sed 's/ -o / /')
echo "Final command to run:"
echo -e "eval "$FIND_CMD -print0" | xargs -0 grep -Hn "^#$TAG" \n"

echo "Instances of $TAG"
eval "$FIND_CMD -print0 | xargs -0 grep -Hn "^#$TAG""
#TODO Create GREP_CMD for colorful output. See Complete_Task: Bash automation Find Tags 
#find "$DIRECTORY" -type f \( -name "*.sh" -o -name "*.R" \) -print0 | xargs -0 grep -Hn "^#$TAG"

NUMBER_OF_INSTANCES=$(eval "$FIND_CMD -print0 | xargs -0 grep -Hn "^#$TAG"" | wc -l )

echo -e "\n$NUMBER_OF_INSTANCES instances of $TAG found"


