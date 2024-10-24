#!/bin/bash

source "../config/slurm_config.sh"

function calculate_indices() {
    local task_id="$1"
    local fastq_count="$2"
    
    local genome_index=$(( (task_id - 1) / fastq_count ))
    local fastq_index=$(( (task_id - 1) % fastq_count ))
    
    echo "${genome_index}:${fastq_index}"
}

function get_output_names() {
    local fastq_path="$1"
    local genome_path="$2"
    
    local fastq_id=$(basename "${fastq_path%.fastq}")
    local genome_name=$(basename "$genome_path" | cut -d_ -f1)
    
    echo "${fastq_id}:${genome_name}"
}

function perform_alignment() {
    local genome_path="$1"
    local fastq_path="$2"
    local output_dir="$3"
    local threads="$4"
    
    local index_base="${genome_path%_refgenome.fna}${FILE_PATTERNS[INDEX_SUFFIX]}"
    local names=$(get_output_names "$fastq_path" "$genome_path")
    local fastq_id=${names%:*}
    local genome_name=${names#*:}
    local output_bam="${output_dir}/${fastq_id}_${genome_name}.bam"
    
    log_info "Starting alignment: $fastq_id to $genome_name"
    
    # Alignment pipeline
    if ! bowtie2 -x "$index_base" \
                 -U "$fastq_path" \
                 -p "$threads" \
                 ${ALIGNMENT_CONFIG[BOWTIE_PARAMS]} |
         samtools view -@ "$threads" -b - |
         samtools sort -@ "$threads" -o "$output_bam" -; then
        log_error "Alignment failed for: $fastq_id to $genome_name"
        return 1
    }
    
    if ! samtools index "$output_bam"; then
        log_error "Index creation failed for: $output_bam"
        return 1
    }
    
    log_info "Completed alignment: $output_bam"
    return 0
}
