#!/bin/bash

source "../config/documentation_config.sh"

function escape_xml_content() {
    local content="$1"
    
    # Remove problematic characters
    content=$(echo "$content" | tr -d '\023')
    
    # Replace XML special characters
    for char in "${!XML_ENTITIES[@]}"; do
        content=${content//$char/${XML_ENTITIES[$char]}}
    done
    
    echo "$content"
}

function write_xml_element() {
    local element="$1"
    local content="$2"
    local indent="${3:-2}"
    local spaces=$(printf '%*s' "$indent" '')
    
    echo "${spaces}<${element}>${content}</${element}>"
}

function create_xml_header() {
    local output_file="$1"
    
    {
        echo '<?xml version="1.0" encoding="UTF-8"?>'
        echo '<repository>'
    } > "$output_file"
}

function close_xml_document() {
    local output_file="$1"
    echo '</repository>' >> "$output_file"
}
