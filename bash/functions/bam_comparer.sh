#!/bin/bash

function generate_output_name() {
    local sample_bam="$1"
    local input_bam="$2"
    local time_id="$3"
    local output_dir="$4"
    
    local sample_name=$(basename "${sample_bam%.bam}" | awk -F'_' '{print $1}')
    local input_name=$(basename "${input_bam%.bam}" | awk -F'_' '{print $1}')
    
    echo "${output_dir}/${time_id}_${sample_name}_${input_name}_bamcomp.bw"
}

function build_bamcompare_command() {
    local sample_bam="$1"
    local input_bam="$2"
    local output_file="$3"
    local threads="$4"
    
    local cmd="bamCompare"
    cmd+=" -b1 \"$sample_bam\""
    cmd+=" -b2 \"$input_bam\""
    cmd+=" -o \"$output_file\""
    cmd+=" --binSize ${BAMCOMPARE_PARAMS[BIN_SIZE]}"
    cmd+=" --normalizeUsing ${BAMCOMPARE_PARAMS[NORMALIZATION]}"
    cmd+=" --scaleFactorsMethod ${BAMCOMPARE_PARAMS[SCALE_METHOD]}"
    cmd+=" --effectiveGenomeSize ${BAMCOMPARE_PARAMS[GENOME_SIZE]}"
    cmd+=" --minMappingQuality ${BAMCOMPARE_PARAMS[MIN_MAPPING_QUALITY]}"
    cmd+=" --operation ${BAMCOMPARE_PARAMS[OPERATION]}"
    cmd+=" --numberOfProcessors $threads"
    
    if [ "${BAMCOMPARE_PARAMS[IGNORE_DUPLICATES]:-true}" = true ]; then
        cmd+=" --ignoreDuplicates"
    fi
    
    for region in "${BAMCOMPARE_PARAMS[IGNORE_REGIONS][@]}"; do
        cmd+=" --ignoreForNormalization $region"
    done
    
    echo "$cmd"
}

function run_comparison() {
    local sample_bam="$1"
    local input_bam="$2"
    local output_file="$3"
    local threads="$4"
    
    log_info "Running BAM comparison"
    log_info "Sample: $sample_bam"
    log_info "Input: $input_bam"
    log_info "Output: $output_file"
    
    local start_time=$(date +%s)
    
    local cmd=$(build_bamcompare_command "$sample_bam" "$input_bam" "$output_file" "$threads")
    
    if ! eval "$cmd"; then
        log_error "Comparison failed"
        return 1
    }
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    
    log_info "Comparison completed in $duration seconds"
    return 0
}
