#!/bin/bash

# Use find to locate .sh and .R files recursively
# Use grep to search for #TODO within the found files
# Handle spaces in filenames and ensure efficient processing

find . -type f \( -name "*.sh" -o -name "*.R" \) -print0 | xargs -0 grep -Hn "^#$1"
