#!/bin/bash

count_string() {
    local search_string="$1"
    local search_dir="${2:-.}"
    local exclude_dirs=(".git" "node_modules" "build" "dist" "renv" ".venv")
    local exclude_files=("*.log" "*.tmp" "*.bak" "*.swp" ".gitignore" ".Rprofile")
    
    find "$search_dir" -type f \
        $(printf -- "-not -path '*%s*' " "${exclude_dirs[@]}") \
        $(printf -- "-not -name '%s' " "${exclude_files[@]}") \
        -exec grep -l "$search_string" {} \; | wc -l
}
