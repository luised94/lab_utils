#!/bin/bash

function find_reference_genomes() {
    local base_dir="${GENOME_CONFIG[BASE_DIR]}"
    local pattern="${GENOME_CONFIG[GENOME_PATTERN]}"
    
    log_info "Searching for reference genomes in: $base_dir"
    
    if [ ! -d "$base_dir" ]; then
        log_error "Reference genome directory not found: $base_dir"
        return 1
    }
    
    find "$base_dir" -type f -name "$pattern"
}

function build_genome_index() {
    local genome_path="$1"
    local index_suffix="${GENOME_CONFIG[INDEX_SUFFIX]}"
    
    if [ ! -f "$genome_path" ]; then
        log_error "Genome file not found: $genome_path"
        return 1
    }
    
    local index_base="${genome_path%_refgenome.fna}${index_suffix}"
    
    log_info "Building index for: $genome_path"
    log_info "Output base: $index_base"
    
    if ! bowtie2-build "$genome_path" "$index_base"; then
        log_error "Index building failed for: $genome_path"
        return 1
    }
    
    log_info "Successfully built index for: $genome_path"
    return 0
}
