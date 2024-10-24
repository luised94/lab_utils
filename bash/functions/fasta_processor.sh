#!/bin/bash

source "../config/genome_config.sh"

function validate_fasta() {
    local file="$1"
    
    log_info "Validating FASTA file: $file"
    
    if [ ! -f "$file" ]; then
        log_error "File not found: $file"
        return 1
    }
    
    if ! grep -q '^>' "$file"; then
        log_error "Invalid FASTA format: No headers found"
        return 1
    }
    
    return 0
}

function create_backup() {
    local source_file="$1"
    local backup_suffix="${2:-${GENOME_FILES[BACKUP_SUFFIX]}}"
    
    local backup_file="${source_file%_refgenome.fna}${backup_suffix}"
    
    log_info "Creating backup: $backup_file"
    
    if [ -f "$backup_file" ]; then
        log_warning "Backup file already exists: $backup_file"
        return 0
    }
    
    if ! cp "$source_file" "$backup_file"; then
        log_error "Failed to create backup"
        return 1
    }
    
    return 0
}

function reformat_headers_awk() {
    local input_file="$1"
    local output_file="$2"
    local old_prefix="${HEADER_PATTERNS[OLD_PREFIX]}"
    local new_prefix="${HEADER_PATTERNS[NEW_PREFIX]}"
    local delimiter="${HEADER_PATTERNS[DELIMITER]}"
    
    log_info "Reformatting headers with AWK"
    
    awk -v old="$old_prefix" -v new="$new_prefix" -v delim="$delimiter" '
        /^>/ {
            gsub(old, new, $6)
            printf(">%s%s\n", $6, $7)
            next
        }
        {
            print $0
        }' "$input_file" | 
    sed "s/${delimiter}.*//" > "$output_file"
}

function reformat_headers_while() {
    local input_file="$1"
    local output_file="$2"
    local old_prefix="${HEADER_PATTERNS[OLD_PREFIX]}"
    local new_prefix="${HEADER_PATTERNS[NEW_PREFIX]}"
    local delimiter="${HEADER_PATTERNS[DELIMITER]}"
    
    log_info "Reformatting headers with while loop"
    
    while IFS= read -r line; do
        if [[ $line == ">"* ]]; then
            echo "${line/$old_prefix/$new_prefix}" | cut -d"$delimiter" -f1
        else
            echo "$line"
        fi
    done < "$input_file" > "$output_file"
}

function verify_conversion() {
    local original="$1"
    local converted="$2"
    
    log_info "Verifying conversion"
    
    local orig_count=$(grep -c '^>' "$original")
    local conv_count=$(grep -c '^>' "$converted")
    
    if [ "$orig_count" -ne "$conv_count" ]; then
        log_error "Header count mismatch: Original=$orig_count, Converted=$conv_count"
        return 1
    }
    
    if grep -q "${HEADER_PATTERNS[OLD_PREFIX]}" "$converted"; then
        log_error "Old prefix still present in converted file"
        return 1
    }
    
    log_info "Conversion verified successfully"
    return 0
}
