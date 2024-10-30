
#!/bin/bash

source "../config/fastq_processing_config.sh"

function setup_output_directories() {
    local base_dir="$1"
    
    log_info "Setting up output directories"
    
    for dir in "${OUTPUT_DIRS[@]}"; do
        local full_path="$base_dir/$dir"
        mkdir -p "$full_path" || {
            log_error "Failed to create directory: $full_path"
            return 1
        }
    done
    
    return 0
}

function find_input_files() {
    local base_dir="$1"
    
    log_info "Finding input FASTQ files"
    
    local exclude_pattern=""
    for pattern in "${FILE_PATTERNS[EXCLUDE_PATTERNS][@]}"; do
        exclude_pattern+=" ! -name \"$pattern\""
    done
    
    eval "find \"$base_dir\" -type f -name \"${FILE_PATTERNS[INPUT]}\" $exclude_pattern"
}

function get_output_paths() {
    local input_file="$1"
    local base_dir="$2"
    
    local filename=$(basename "$input_file")
    local processed_name="${FILE_PATTERNS[OUTPUT_PREFIX]}${filename}"
    
    echo "${base_dir}/${OUTPUT_DIRS[PROCESSED]}/${processed_name}:${base_dir}/${OUTPUT_DIRS[LOGS]}/${processed_name%.fastq}.json"
}

function build_fastp_command() {
    local input_file="$1"
    local output_file="$2"
    local json_file="$3"
    
    local params=("${FASTP_PARAMS[BASE_PARAMS][@]}")
    
    if [[ $input_file =~ "Eaton" ]]; then
        log_info "Using Eaton-specific parameters"
        params+=("${FASTP_PARAMS[EATON][@]}")
    else
        params+=("${FASTP_PARAMS[STANDARD][@]}")
    fi
    
    echo "fastp -i \"$input_file\" -o \"$output_file\" --json \"$json_file\" ${params[*]}"
}

function process_fastq_file() {
    local input_file="$1"
    local base_dir="$2"
    
    log_info "Processing file: $input_file"
    
    local paths=$(get_output_paths "$input_file" "$base_dir")
    local output_file=${paths%:*}
    local json_file=${paths#*:}
    
    local cmd=$(build_fastp_command "$input_file" "$output_file" "$json_file")
    
    log_info "Executing: $cmd"
    if ! eval "$cmd"; then
        log_error "Processing failed for: $input_file"
        return 1
    }
    
    log_info "Successfully processed: $input_file"
    return 0
}
