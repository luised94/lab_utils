#!/bin/bash

source "../functions/xml_handler.sh"

function get_repository_root() {
    if command -v git &>/dev/null && git rev-parse --is-inside-work-tree &>/dev/null; then
        git rev-parse --show-toplevel
        return 0
    fi
    
    log_warning "Git not available or not in a repository"
    echo "."
    return 1
}

function build_exclusion_pattern() {
    local patterns=()
    for pattern in "${EXCLUDED_PATTERNS[@]}"; do
        patterns+=("-not" "-path" "*/$pattern/*")
    done
    echo "${patterns[@]}"
}

function process_file() {
    local file="$1"
    local output_file="$2"
    
    local file_name=$(basename "$file")
    local file_path=$(dirname "$file")
    local mime_type=$(file -b --mime-type "$file")
    
    {
        echo "<file>"
        write_xml_element "name" "$(escape_xml_content "$file_name")"
        write_xml_element "path" "$(escape_xml_content "$file_path")"
        write_xml_element "type" "$mime_type"
        
        echo "  <content><![CDATA["
        if [ -r "$file" ]; then
            cat "$file" || log_error "Failed to read: $file"
        else
            log_error "File not readable: $file"
            echo "Error: Unable to read $file"
        fi
        echo "]]></content>"
        echo "</file>"
    } >> "$output_file"
}

function process_file_type() {
    local file_type="$1"
    local output_file="$2"
    local root_dir="$3"
    
    log_info "Processing .$file_type files"
    
    local exclusions=$(build_exclusion_pattern)
    local find_cmd="find \"$root_dir\" -type f -name \"*.$file_type\" $exclusions -print0"
    
    while IFS= read -r -d '' file; do
        log_info "Processing: $file"
        process_file "$file" "$output_file"
    done < <(eval "$find_cmd")
}
