#STATUS:
#!/bin/bash

find . -maxdepth 1 -type f -name "slurm*out" -exec echo {} \;

# Define the directory to search and the pattern to match
SEARCH_DIR="/path/to/search"
PATTERN="*.tmp"  # Example pattern, adjust as needed

# Perform the dry run: list files that would be deleted
echo "Performing dry run..."
find "$SEARCH_DIR" -type f -name "$PATTERN" -exec echo "Would delete: {}" \;

# Ask the user for confirmation to proceed with actual deletion
read -p "Proceed with deleting these files? (y/n): " confirm && [[ $confirm == [yY] ]] || exit 1

# If confirmed, proceed with deletion
echo "Deleting files..."
find "$SEARCH_DIR" -type f -name "$PATTERN" -exec rm {} \;

echo "Files have been deleted.
