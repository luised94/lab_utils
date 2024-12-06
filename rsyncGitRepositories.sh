#STATUS: KEEP.
# Move to code_review_holder.
#!/bin/bash
#if [ $# -ne 1 ]; then
echo "Usage: $0 "
#fi
# Set default values for source and destination directories
SOURCE_DIR="$HOME"
DEST_DIR="/mnt/c/Users/${WINDOWS_USER}/Dropbox (MIT)/"

# Check if the source and destination directories exist
if [ ! -d "$SOURCE_DIR" ]; then
    echo "Error: Source directory '$SOURCE_DIR' does not exist."
    exit 1
fi

if [ ! -d "$DEST_DIR" ]; then
    printf "Error: Destination directory '%s' does not exist.\n" "$DEST_DIR"
    printf "Please ensure that the Dropbox folder is mounted and accessible.\n"
    exit 1
fi

# Print source and destination directories for debugging
printf "Source directory: %s\n" "$SOURCE_DIR"
printf "Destination directory: %s\n" "$DEST_DIR"

find "$SOURCE_DIR" -maxdepth 1 -type d -print0 | while IFS= read -r -d $'\0' dir; do

  # Check if the directory is a Git repository
  if [ -d "$dir/.git" -a "$dir" != "$HOME/.nvm" -a "$dir" != "$HOME/cytolib" ]; then
    repo_name=$(basename "$dir")

    printf "Found Git repository: %s\n" "$dir"
    printf "Syncing to: %s\n" "$DEST_DIR"

    # Rsync the Git repository, excluding the .git directory
    rsync -av --delete --exclude=.git/ "$dir" "$DEST_DIR"

    if [ $? -eq 0 ]; then
      printf "Successfully synced %s.\n" "$repo_name"
    else
      printf "Error syncing %s. Please check the logs.\n" "$repo_name"
    fi

  fi

done
