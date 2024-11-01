#!/bin/bash

source "../config/genome_config.sh"

function extract_organism_name() {
    local report_file="$1"
    
    if [ ! -f "$report_file" ]; then
        log_error "Assembly report not found: $report_file"
        return 1
    fi
    
    log_info "Extracting organism name from: $report_file"
    
    local organism_name=$(grep -m 1 -o '"organismName":"[^"]*' "$report_file" | 
                         cut -d '"' -f4 | 
                         sed 's/ //g')
    
    if [ -z "$organism_name" ]; then
        log_error "Failed to extract organism name"
        return 1
    fi
    
    echo "$organism_name"
}

function reorganize_genome_files() {
    local source_dir="$1"
    local organism_name="$2"
    
    log_info "Reorganizing files for: $organism_name"
    
    # Move all files to root directory
    find "$source_dir/" -type f -exec mv -v {} "$source_dir/" \; || {
        log_error "Failed to move files in: $source_dir"
        return 1
    }
    
    # Rename CDS file
    local cds_source="${source_dir}/${GENOME_NAMING[ORIGINAL_CDS]}"
    local cds_target="${source_dir}/${GENOME_NAMING[CDS_FILE]}"
    if [ -f "$cds_source" ]; then
        mv -v "$cds_source" "$cds_target" || log_warning "Failed to rename CDS file"
    fi
    
    # Rename genomic file
    find "$source_dir" -name "${GENOME_NAMING[GENOMIC_PATTERN]}" \
        -exec mv -v {} "$source_dir/${organism_name}${GENOME_NAMING[REFGENOME_SUFFIX]}" \;
    
    # Cleanup
    rm -rf "$source_dir/${GENOME_PATHS[NCBI_DATA]}"
    
    return 0
}

function standardize_chromosome_names() {
    local genome_file="$1"
    local output_file="$2"
    
    log_info "Standardizing chromosome names in: $genome_file"
    
    awk -v old="${CHROMOSOME_NAMING[PATTERN_OLD]}" \
        -v new="${CHROMOSOME_NAMING[PATTERN_NEW]}" \
        -v sep="${CHROMOSOME_NAMING[SEPARATOR]}" \
        '/^>/ {
            gsub(old, new, $6)
            printf(">%s%s\n", $6, $7)
            next
        }
        !/^>/ {
            print $0
        }' "$genome_file" | 
    sed "s/${CHROMOSOME_NAMING[SEPARATOR]}.*//" > "$output_file"
}
