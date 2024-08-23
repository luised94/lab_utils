#!/bin/bash
if [ $# -gt 2 ]; then
  echo "Usage: $0 [<source_directory>] [<destination_directory>]"
  exit 1
fi
# Set default values for source and destination directories
SOURCE_DIR="$HOME"
DEST_DIR="/mnt/c/Users/${1}/Dropbox (MIT)/"

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


