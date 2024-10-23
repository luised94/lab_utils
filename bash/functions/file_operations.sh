#!/bin/bash
count_string() {
    local search_string="$1"
    local search_dir="${2:-.}"
    local exclude_dirs=(".git" "node_modules" "build" "dist" "renv" ".venv")
    local exclude_files=("*.log" "*.tmp" "*.bak" "*.swp" ".gitignore" ".Rprofile")
    
    # Build exclude arguments properly
    local exclude_args=()
    for dir in "${exclude_dirs[@]}"; do
        exclude_args+=(-not -path "*/${dir}/*")
    done
    for file in "${exclude_files[@]}"; do
        exclude_args+=(-not -name "${file}")
    done

    echo "Searching for files containing: '$search_string'"
    echo "----------------------------------------"
    
    # Files containing the string
    echo "Files containing the string:"
    local files_with=$(find "$search_dir" -type f "${exclude_args[@]}" -exec grep -l "$search_string" {} \;)
    local count_with=$(echo "$files_with" | grep -c "^" || echo 0)
    
    if [ $count_with -gt 0 ]; then
        echo "$files_with" | sed 's/^/  /'
    else
        echo "  None found"
    fi
    
    echo -e "\nFiles missing the string:"
    local files_without=$(find "$search_dir" -type f "${exclude_args[@]}" -exec grep -L "$search_string" {} \;)
    local count_without=$(echo "$files_without" | grep -c "^" || echo 0)
    
    if [ $count_without -gt 0 ]; then
        echo "$files_without" | sed 's/^/  /'
    else
        echo "  None found"
    fi
    
    echo -e "\nSummary:"
    echo "  Files containing string: $count_with"
    echo "  Files missing string: $count_without"
    echo "  Total files checked: $((count_with + count_without))"
}
