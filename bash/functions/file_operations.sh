#!/bin/bash
count_string() {
    # Usage function
    usage() {
        cat << EOF
Usage: count_string [OPTIONS] <search_string> [directory]
Search for string occurrences in files with detailed reporting.

Options:
    -h, --help                 Show this help message
    -e, --exclude-dir DIR      Additional directory to exclude (can be used multiple times)
    -f, --exclude-file FILE    Additional file pattern to exclude (can be used multiple times)
    -v, --verbose             Enable verbose output
    -q, --quiet               Suppress all output except final counts
    --no-default-excludes     Don't use default exclusion patterns
    --max-depth N             Maximum directory depth to search

Examples:
    count_string "TODO" ./src
    count_string -e "tests" -e "docs" "FIXME" .
    count_string -q "deprecated" ./project
EOF
    }

    # Default configuration
    local default_exclude_dirs=(".git" "node_modules" "build" "dist" "renv" ".venv")
    local default_exclude_files=("*.md" "*.txt" "*init.sh" "*renv.lock" "*.log" "*.tmp" "*.bak" "*.swp" "*.gitignore" "*.Rprofile")
    local additional_exclude_dirs=()
    local additional_exclude_files=()
    local verbose=0
    local quiet=0
    local use_default_excludes=1
    local max_depth=""
    local search_string=""
    local search_dir="."

    # Parse options
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                usage
                return 0
                ;;
            -e|--exclude-dir)
                if [[ -z "$2" ]]; then
                    echo "Error: --exclude-dir requires a directory argument" >&2
                    return 1
                fi
                additional_exclude_dirs+=("$2")
                shift 2
                ;;
            -f|--exclude-file)
                if [[ -z "$2" ]]; then
                    echo "Error: --exclude-file requires a file pattern argument" >&2
                    return 1
                fi
                additional_exclude_files+=("$2")
                shift 2
                ;;
            -v|--verbose)
                verbose=1
                shift
                ;;
            -q|--quiet)
                quiet=1
                shift
                ;;
            --no-default-excludes)
                use_default_excludes=0
                shift
                ;;
            --max-depth)
                if [[ -z "$2" ]] || ! [[ "$2" =~ ^[0-9]+$ ]]; then
                    echo "Error: --max-depth requires a numeric argument" >&2
                    return 1
                fi
                max_depth="-maxdepth $2"
                shift 2
                ;;
            -*)
                echo "Error: Unknown option: $1" >&2
                usage
                return 1
                ;;
            *)
                if [[ -z "$search_string" ]]; then
                    search_string="$1"
                else
                    search_dir="$1"
                fi
                shift
                ;;
        esac
    done

    # Validate required arguments
    if [[ -z "$search_string" ]]; then
        echo "Error: Search string is required" >&2
        usage
        return 1
    fi

    # Validate directory
    if [[ ! -d "$search_dir" ]]; then
        echo "Error: Directory '$search_dir' does not exist" >&2
        return 1
    fi

    # Build exclude arguments
    local exclude_args=()
    
    if ((use_default_excludes)); then
        for dir in "${default_exclude_dirs[@]}"; do
            exclude_args+=(-not -path "*/${dir}/*")
        done
        for file in "${default_exclude_files[@]}"; do
            exclude_args+=(-not -name "${file}")
        done
    fi

    for dir in "${additional_exclude_dirs[@]}"; do
        exclude_args+=(-not -path "*/${dir}/*")
    done
    for file in "${additional_exclude_files[@]}"; do
        exclude_args+=(-not -name "${file}")
    done

    # Temporary files for results
    local tmp_dir=$(mktemp -d)
    local files_with="$tmp_dir/with.txt"
    local files_without="$tmp_dir/without.txt"
    trap 'rm -rf "$tmp_dir"' EXIT

    # Execute find command with proper error handling
    if ((verbose)); then
        echo "Executing find command..."
        echo "find $search_dir $max_depth -type f ${exclude_args[@]}"
    fi

    # Find and categorize files
    find "$search_dir" $max_depth -type f "${exclude_args[@]}" -print0 2>/dev/null | \
        while IFS= read -r -d $'\0' file; do
            if grep -q "$search_string" "$file" 2>/dev/null; then
                echo "$file" >> "$files_with"
            else
                echo "$file" >> "$files_without"
            fi
        done

    # Count results
    local count_with=$(wc -l < "$files_with" || echo 0)
    local count_without=$(wc -l < "$files_without" || echo 0)
    local total=$((count_with + count_without))

    # Output results
    if ((! quiet)); then
        echo "Searching for: '$search_string' in $search_dir"
        echo "----------------------------------------"
        
        echo -e "\nFiles containing the string:"
        if [[ -s "$files_with" ]]; then
            sed 's/^/  /' "$files_with"
        else
            echo "  None found"
        fi
        
        echo -e "\nFiles missing the string:"
        if [[ -s "$files_without" ]]; then
            sed 's/^/  /' "$files_without"
        else
            echo "  None found"
        fi
        
        echo -e "\nSummary:"
        echo "  Files containing string: $count_with"
        echo "  Files missing string: $count_without"
        echo "  Total files checked: $total"
    else
        echo "$count_with"
    fi

    return 0
}
