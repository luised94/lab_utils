
determine_next_experiment_id() {
    #Find all experiment*.md files, extract second element using _ delimiter. If it is three digit number,
    # select the maximum. if none is found, outputs 001.
    echo $(find . -maxdepth 1 -type f -name "*.sh" | sed 's/\.\///g' | awk -F'_' '
        $1 ~ /^[0-9]{3}$/ {
            if ($1 > max) max = $1
        }
        END {
        printf "%03d", (max == "" ? 1 : max + 1)
        }
    ')
}

# Create the name of the file, find the template location (./templates/ relative to the script location)
# Used sed to replace the tags on the template file with appropriate values. 
#TODO Have to add descriptive name section that reads in name.
#TODO HAve to figure out how to deal with connected files, files with multiple stages, and use template for more specific experiments.
create_new_experiment() {
    local date_of_creation=$(date +%Y%m%d)
    experiment_index=$(determine_next_experiment_id)
    local filename="${date_of_creation}_${experiment_index}_experiment.md"
    SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

    #template_file="${SCRIPT_DIR}/templates/experiment_template.md" 
    #sed "s/{{EXPERIMENT_NAME}}/trial/g; s/{{DATE}}/${date_of_creation}/g" "$template_file" > "$filename"
    #nvim $filename
    echo $experiment_index
}

#!/bin/bash

source "../config/sra_config.sh"

function validate_input() {
    local download_dir="$1"
    
    log_info "Validating input parameters"
    
    if [ -z "$download_dir" ]; then
        log_error "Download directory not specified"
        return 1
    fi
    
    local full_path="${SRA_CONFIG[DATA_DIR]}/$download_dir"
    
    if [ ! -d "$full_path" ]; then
        log_info "Creating directory: $full_path"
        mkdir -p "$full_path" || {
            log_error "Failed to create directory: $full_path"
            return 1
        }
    fi
    
    echo "$full_path"
}

function construct_download_url() {
    local accession="$1"
    local base_url="${SRA_CONFIG[BASE_URL]}"
    
    echo "${base_url}${accession:0:6}/${accession}/${accession}.fastq.gz"
}

function verify_url() {
    local url="$1"
    
    log_info "Verifying URL: $url"
    
    if ! curl --head --silent --fail "$url" >/dev/null; then
        log_error "URL not accessible: $url"
        return 1
    fi
    
    return 0
}

function download_file() {
    local url="$1"
    local output_file="$2"
    
    log_info "Downloading: $url -> $output_file"
    
    if ! wget --quiet --show-progress --output-document="$output_file" "$url"; then
        log_error "Download failed: $url"
        return 1
    fi
    
    log_info "Download complete: $output_file"
    return 0
}

function concatenate_files() {
    local output_dir="$1"
    local output_file="$2"
    local files=("${@:3}")
    
    log_info "Concatenating files to: $output_file"
    
    for file in "${files[@]}"; do
        if [ ! -f "$output_dir/$file" ]; then
            log_error "File not found: $file"
            return 1
        fi
        cat "$output_dir/$file" >> "$output_file" || {
            log_error "Failed to concatenate: $file"
            return 1
        }
    done
    
    log_info "Concatenation complete"
    return 0
}

function decompress_file() {
    local file="$1"
    
    log_info "Decompressing: $file"
    
    if ! gunzip "$file"; then
        log_error "Decompression failed: $file"
        return 1
    fi
    
    log_info "Decompression complete"
    return 0
}

#!/bin/bash

source "../config/bam_comparison_config.sh"

function get_bam_pairs() {
    local base_dir="$1"
    local task_id="$2"
    local r_script="${R_CONFIG[SCRIPT_PATH]}"
    
    log_info "Getting BAM pairs from R script"
    
    if [ ! -f "$r_script" ]; then
        log_error "R script not found: $r_script"
        return 1
    fi
    
    local pairs=$(Rscript "$r_script" "$base_dir" "$task_id" | grep -v '^NULL$')
    
    if [ -z "$pairs" ]; then
        log_error "No output from R script"
        return 1
    fi
    
    echo "$pairs"
}

function validate_bam_pairs() {
    local -a pairs=("$@")
    
    log_info "Validating BAM pairs"
    
    if [ ${#pairs[@]} -ne 2 ]; then
        log_error "Expected 2 BAM files, got ${#pairs[@]}"
        return 1
    fi
    
    for bam in "${pairs[@]}"; do
        if [ ! -f "$bam" ]; then
            log_error "BAM file not found: $bam"
            return 1
        fi
    done
    
    return 0
}

#!/bin/bash

source "../functions/ngs_file_manager.sh"

function setup_qc_directory() {
    local base_dir="$1"
    local qc_dir="${base_dir}/${QC_CONFIG[OUTPUT_DIR]}"
    
    log_info "Setting up QC directory: $qc_dir"
    mkdir -p "$qc_dir" || {
        log_error "Failed to create QC directory: $qc_dir"
        return 1
    }
    
    echo "$qc_dir"
}

function run_fastqc() {
    local file_path="$1"
    local output_dir="$2"
    
    if [ ! -f "$file_path" ]; then
        log_error "File not found: $file_path"
        return 1
    fi
    
    log_info "Running FASTQC on: $file_path"
    if ! fastqc --outdir="$output_dir" "$file_path"; then
        log_error "FASTQC failed for: $file_path"
        return 1
    fi
    
    log_info "FASTQC completed for: $file_path"
    return 0
}

#!/bin/bash

source "../config/file_management_config.sh"

function build_file_patterns() {
    local patterns=()
    for type in "${!FILE_TYPES[@]}"; do
        for ext in ${FILE_TYPES[$type]}; do
            patterns+=("-o" "-name" "*.$ext")
        done
    done
    echo "${patterns[@]:1}" # Remove first -o
}

function check_disk_space() {
    local target_dir="$1"
    local min_space="${2:-${DEFAULTS[MIN_FREE_SPACE_GB]}}"
    
    local available_space=$(df -BG "$target_dir" | awk 'NR==2 {print $4}' | sed 's/G//')
    if (( available_space < min_space )); then
        log_error "Insufficient disk space. Available: ${available_space}GB, Required: ${min_space}GB"
        return 1
    if
    return 0
}

function validate_directories() {
    local target_dir="$1"
    local search_dir="$2"
    
    if [[ ! -d "$search_dir" ]]; then
        log_error "Search directory does not exist: $search_dir"
        return 1
    fi
    
    if [[ ! -d "$target_dir" ]]; then
        log_info "Creating target directory: $target_dir"
        mkdir -p "$target_dir"
    fi
    
    if [[ ! -w "$target_dir" ]]; then
        log_error "Target directory not writable: $target_dir"
        return 1
    fi
}

function move_ngs_files() {
    local target_dir="$1"
    local search_dir="${2:-.}"
    local max_depth="${3:-${DEFAULTS[MAX_DEPTH]}}"
    local batch_size="${4:-${DEFAULTS[BATCH_SIZE]}}"
    
    log_info "Starting NGS file movement operation"
    
    validate_directories "$target_dir" "$search_dir" || return 1
    check_disk_space "$target_dir" || return 1
    
    local file_patterns=($(build_file_patterns))
    
    find "$search_dir" \
        -maxdepth "$max_depth" \
        -type f \
        \( "${file_patterns[@]}" \) \
        \( -path "*/code/*" -o -path "*/script*/*" \) \
        -print0 | 
    while IFS= read -r -d '' file; do
        local file_type=$(determine_file_type "$file")
        local target_subdir="$target_dir/$file_type"
        
        mkdir -p "$target_subdir"
        log_info "Moving $file to $target_subdir"
        mv "$file" "$target_subdir/" || log_error "Failed to move: $file"
    done
}

function determine_file_type() {
    local file="$1"
    local extension="${file##*.}"
    
    for type in "${!FILE_TYPES[@]}"; do
        if [[ " ${FILE_TYPES[$type]} " =~ " $extension " ]]; then
            echo "$type"
            return
        }
    done
    echo "OTHER"
}

function analyze_file_distribution() {
    local search_dir="${1:-.}"
    
    log_info "Analyzing file distribution in $search_dir"
    
    find "$search_dir" \
        -type d \( -path '*/lib/*' -o -path '*/renv/*' -o -path '*/git/*' \) -prune \
        -o -type f -exec du -a {} + | 
        sort -nr | 
        head -n 1000 | 
        awk -F/ '{print $NF}' | 
        rev | cut -d. -f1 | rev | 
        sort | uniq -c | 
        sort -nr
}

# Add to existing file management functions
function find_fastq_files() {
    local base_dir="$1"
    local exclude_pattern="${QC_CONFIG[EXCLUDE_PATTERNS]}"
    
    log_info "Finding FASTQ files in: $base_dir"
    find "$base_dir" -type f -name "*.fastq" ! \( -name "*unmapped*" -o -name "processed_*" \)
}

function find_processed_fastq() {
    local base_dir="$1"
    local pattern="${QC_CONFIG[PROCESSED_PATTERN]}"
    
    log_info "Finding processed FASTQ files in: $base_dir"
    find "$base_dir" -type f -name "$pattern"
}

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

#!/bin/bash

source "../config/environment_config.sh"

function validate_directory() {
    local dir="$1"
    
    log_info "Validating directory: $dir"
    
    if [ ! -d "$dir" ]; then
        log_error "Directory not found: $dir"
        return 1
    fi
    
    if [ ! -r "$dir" ]; then
        log_error "Directory not readable: $dir"
        return 1
    fi
    
    return 0
}

function find_function_files() {
    local base_dir="$1"
    local pattern="${FILE_PATTERNS[FUNCTIONS]}"
    
    log_info "Finding function files"
    
    local exclude_pattern=""
    for pattern in "${FILE_PATTERNS[EXCLUDE_PATTERNS][@]}"; do
        exclude_pattern+=" ! -name \"$pattern\""
    done
    
    eval "find \"$base_dir\" -type f -name \"$pattern\" $exclude_pattern"
}

function validate_function_file() {
    local file="$1"
    
    log_info "Validating function file: $file"
    
    if [ ! -f "$file" ]; then
        log_error "File not found: $file"
        return 1
    fi
    
    if [ ! -r "$file" ]; then
        log_error "File not readable: $file"
        return 1
    fi
    
    # Optional: Add syntax check
    if command -v bash > /dev/null; then
        if ! bash -n "$file"; then
            log_error "Syntax error in: $file"
            return 1
        fi
    fi
    
    return 0
}

function load_function_file() {
    local file="$1"
    
    log_info "Loading functions from: $file"
    
    if ! source "$file" 2>/dev/null; then
        log_error "Failed to source: $file"
        return 1
    fi
    
    return 0
}

function load_priority_functions() {
    local base_dir="$1"
    
    log_info "Loading priority functions"
    
    for file in "${LOAD_ORDER[PRIORITY][@]}"; do
        local full_path="$base_dir/$file"
        if [ -f "$full_path" ]; then
            load_function_file "$full_path" || continue
        else
            log_warning "Priority file not found: $file"
        fi
    done
}

function load_remaining_functions() {
    local base_dir="$1"
    
    log_info "Loading remaining functions"
    
    while IFS= read -r -d $'\0' file; do
        # Skip priority and optional files
        local basename=$(basename "$file")
        if [[ " ${LOAD_ORDER[PRIORITY][@]} ${LOAD_ORDER[OPTIONAL][@]} " =~ " $basename " ]]; then
            continue
        fi
        
        load_function_file "$file" || continue
        
    done < <(find_function_files "$base_dir")
}

#!/bin/bash
# bash/functions/fastq_processor.sh

source "$HOME/lab_utils/bash/functions/logging_utils.sh"

#' Validate FASTQ Processing Input
#' @param experiment_id Character Experiment identifier
#' @param log_file Character Log file path
#' @return String Experiment directory path
validate_fastq_input() {
    local experiment_id="$1"
    local log_file="$2"
    
    log_info "Validating input for experiment: $experiment_id" "$log_file"
    
    # Validate experiment ID format
    if [[ ! "$experiment_id" =~ ^[0-9]{6}Bel$ ]]; then
        log_error "Invalid experiment ID format: $experiment_id" "$log_file"
        return 1
    fi
    
    local experiment_dir="$HOME/data/$experiment_id"
    if [[ ! -d "$experiment_dir" ]]; then
        log_error "Experiment directory not found: $experiment_dir" "$log_file"
        return 1
    fi
    
    echo -n "$experiment_dir"
}

#' Setup FASTQ Processing Directories
#' @param experiment_dir Character Experiment directory
#' @param log_file Character Log file path
#' @return String Output directory path
setup_fastq_directories() {
    local experiment_dir="$1"
    local log_file="$2"
    

    local output_dir="${experiment_dir}/fastq"
    log_info "Setting up output directory: $output_dir" "$log_file"
    
    mkdir -p "$output_dir" || {
        log_error "Failed to create output directory: $output_dir" "$log_file"
        return 1
    }
    
    echo -n "$output_dir"
}

#' Find FASTQ Files
#' @param experiment_dir Character Experiment directory
#' @param log_file Character Log file path
#' @return Array FASTQ file paths

find_fastq_files() {
    local experiment_dir="$1"
    local log_file="$2"
    
    log_info "Searching for FASTQ files in: $experiment_dir/${PROJECT_CONFIG[FASTQ_DIR]}" "$log_file"
    
    # Validate directory exists
    if [[ ! -d "$experiment_dir/${PROJECT_CONFIG[FASTQ_DIR]}" ]]; then
        log_error "FASTQ directory not found: $experiment_dir/${PROJECT_CONFIG[FASTQ_DIR]}" "$log_file"
        return 1
    fi
    
    # Build exclude patterns array
    local -a exclude_patterns=()
    for pattern in ${PROJECT_CONFIG[FASTQ_EXCLUDE]}; do
        exclude_patterns+=( "!" "-name" "$pattern" )
    done
    
    # Execute find with error handling
    local find_output
    if ! find_output=$(find "$experiment_dir/${PROJECT_CONFIG[FASTQ_DIR]}" -type f \
        -name "${PROJECT_CONFIG[FASTQ_PATTERN]}" \
        "${exclude_patterns[@]}" 2>&1 | sort); then
        log_error "Find command failed: $find_output" "$log_file"
        return 1
    fi
    
    # Check if any files were found
    if [[ -z "$find_output" ]]; then
        log_warning "No FASTQ files found matching pattern: ${PROJECT_CONFIG[FASTQ_PATTERN]}" "$log_file"
        return 0
    fi
    
    echo -n "$find_output"
}

#' Process FASTQ Files
#' @param experiment_dir Character Experiment directory
#' @param output_dir Character Output directory
#' @param log_file Character Log file path
#' @return Integer 0 if successful
consolidate_fastq_files_by_id() {
    local experiment_dir="$1"
    local output_dir="$2"
    local log_file="$3"

    local -a fastq_files
    mapfile -t fastq_files < <(find_fastq_files "$experiment_dir" "$log_file")

    local initial_count=${#fastq_files[@]}
    log_info "Found $initial_count FASTQ files to process" "$log_file"

    for file in "${fastq_files[@]}"; do
        local basename=$(basename "$file")
        # Split by both - and _ and look for 5-6 digit pattern
        local id=""
        local parts
        IFS='_-' read -ra parts <<< "$basename"
        for part in "${parts[@]}"; do
            if [[ $part =~ ^[0-9]{5,6}$ ]]; then
                id="$part"
                break
            fi
        done

        if [[ -z "$id" ]]; then
            log_error "Could not extract ID from filename: $basename" "$log_file"
            return 1
        fi
        log_info "Extract id: $id"
        local output_file="$output_dir/${id}${PROJECT_CONFIG[FASTQ_SUFFIX]}"
        log_info "Processing: $basename -> $(basename "$output_file")" "$log_file"

        # Use temporary file for safety
        local temp_output="${output_file}.tmp"
        if ! cat "$file" >> "$temp_output"; then
            log_error "Failed to process file: $file" "$log_file"
            rm -f "$temp_output"  # Clean up temp file
            return 1
        fi
        
        # Move temp file to final location
        if ! mv "$temp_output" "$output_file"; then
            log_error "Failed to finalize output file: $output_file" "$log_file"
            return 1
        fi
        
        # Only remove original after successful processing
        if ! rm "$file"; then
            log_error "Failed to remove original file: $file" "$log_file"
            return 1
        fi
    
        log_info "Successfully processed and removed: $basename" "$log_file"
    done
    
    return 0
}

#!/bin/bash
#' Move fastq files inside directories to current working directory.
#' @param target_dir Character Target directory
#' @param log_file Character Log file path
#' @return Integer 0 if successful, 1 otherwise
move_fastq_files_to_current_directory() {
    local target_dir="$1"
    local log_file="$2"
    log_info "Organizing FASTQ files in: $target_dir" "$log_file"
    # Store current directory
    local current_dir=$(pwd)
    # Change to target directory
    cd "$target_dir" || {
        log_error "Failed to access directory: $target_dir" "$log_file"
        return 1
    }
    # Move files to root of target directory
    find . -type f -name "*.fastq" -exec mv {} . \; || {
        log_error "Failed to move FASTQ files" "$log_file"
        cd "$current_dir"
        return 1
    }
    # Return to original directory
    cd "$current_dir"
    log_info "FASTQ files organized successfully" "$log_file"
    return 0
}

clean_experiment_directory() {
    local target_dir="$1"
    local log_file="$2"
    log_info "Starting cleanup process in: $target_dir" "$log_file"
    # Store current directory
    local current_dir=$(pwd)
    # Verify we're in the correct directory structure
    if [[ "$target_dir" != "${BMC_CONFIG[TARGET_FS]}"* ]]; then
        log_error "Invalid target directory: $target_dir" "$log_file"
        log_error "Must be under: ${BMC_CONFIG[TARGET_FS]}" "$log_file"
        return 1
    fi
    # Change to target directory
    cd "$target_dir" || {
        log_error "Failed to access directory: $target_dir" "$log_file"
        return 1
    }

    # Process cleanup patterns
    local -a patterns=(${BMC_CONFIG[CLEANUP_DIRS]} ${BMC_CONFIG[CLEANUP_FILES]})
    for pattern in "${patterns[@]}"; do
        remove_files_safely "$pattern" "$log_file"
    done
    # Return to original directory
    cd "$current_dir"
    
    log_info "Cleanup process completed" "$log_file"
    return 0
}

#!/bin/bash

source "../config/coverage_config.sh"

function setup_coverage_directories() {
    local base_dir="$1"
    
    log_info "Setting up coverage directories"
    
    for dir in "${OUTPUT_DIRS[@]}"; do
        local full_path="$base_dir/$dir"
        mkdir -p "$full_path" || {
            log_error "Failed to create directory: $full_path"
            return 1
        }
    done
    
    return 0
}

function find_bam_files() {
    local base_dir="$1"
    local pattern="${FILE_PATTERNS[BAM_SUFFIX]}"
    
    log_info "Finding BAM files"
    
    local bam_files=$(find "$base_dir" -type f -name "*${pattern}" | sort)
    
    if [ -z "$bam_files" ]; then
        log_error "No BAM files found matching pattern: $pattern"
        return 1
    }
    
    echo "$bam_files"
}

function generate_output_name() {
    local bam_path="$1"
    local time_id="$2"
    local output_dir="$3"
    
    local basename=$(basename "${bam_path%.bam}")
    echo "${output_dir}/${time_id}_${basename}${FILE_PATTERNS[BIGWIG_SUFFIX]}"
}

function build_bamcoverage_command() {
    local input_file="$1"
    local output_file="$2"
    
    local cmd="bamCoverage"
    cmd+=" -b \"$input_file\""
    cmd+=" -o \"$output_file\""
    cmd+=" --binSize ${COVERAGE_PARAMS[BIN_SIZE]}"
    cmd+=" --normalizeUsing ${COVERAGE_PARAMS[NORMALIZATION]}"
    cmd+=" --minMappingQuality ${COVERAGE_PARAMS[MIN_MAPPING_QUALITY]}"
    
    if [ "${COVERAGE_PARAMS[IGNORE_DUPLICATES]}" = true ]; then
        cmd+=" --ignoreDuplicates"
    fi
    
    echo "$cmd"
}

function process_bam_coverage() {
    local bam_file="$1"
    local output_file="$2"
    
    log_info "Processing coverage for: $bam_file"
    log_info "Output: $output_file"
    
    local cmd=$(build_bamcoverage_command "$bam_file" "$output_file")
    
    log_info "Executing: $cmd"
    if ! eval "$cmd"; then
        log_error "Coverage generation failed for: $bam_file"
        return 1
    fi
    
    log_info "Successfully generated coverage for: $bam_file"
    return 0
}

#!/bin/bash

source "../config/quality_control_config.sh"

function setup_qc_directories() {
    local base_dir="$1"
    
    log_info "Setting up QC directories"
    
    for dir in "${QC_DIRS[@]}"; do
        local full_path="$base_dir/$dir"
        mkdir -p "$full_path" || {
            log_error "Failed to create directory: $full_path"
            return 1
        }
    done
    
    return 0
}

function find_bam_files() {
    local base_dir="$1"
    
    log_info "Finding BAM files"
    
    local bam_files=$(find "$base_dir" -type f -name "*.bam")
    
    if [ -z "$bam_files" ]; then
        log_error "No BAM files found in: $base_dir"
        return 1
    }
    
    echo "$bam_files"
}

function get_output_basename() {
    local bam_path="$1"
    
    echo "$(basename "${bam_path%.bam}")"
}

function run_flagstat() {
    local bam_file="$1"
    local output_file="$2"
    
    log_info "Running flagstat on: $bam_file"
    
    if ! samtools flagstat ${SAMTOOLS_PARAMS[FLAGSTAT]} "$bam_file" > "$output_file"; then
        log_error "Flagstat failed for: $bam_file"
        return 1
    }
    
    return 0
}

function run_quickcheck() {
    local bam_file="$1"
    local output_file="$2"
    
    log_info "Running quickcheck on: $bam_file"
    
    if samtools quickcheck "$bam_file"; then
        echo -e 'QUICKCHECK\tTRUE' > "$output_file"
    else
        echo -e 'QUICKCHECK\tFALSE' > "$output_file"
        log_warning "Quickcheck failed for: $bam_file"
    fi
    
    return 0
}

function run_stats() {
    local bam_file="$1"
    local output_file="$2"
    
    log_info "Running stats on: $bam_file"
    
    if ! samtools stats ${SAMTOOLS_PARAMS[STATS]} "$bam_file" > "$output_file"; then
        log_error "Stats failed for: $bam_file"
        return 1
    }
    
    return 0
}

function process_bam_file() {
    local bam_file="$1"
    local qc_dir="$2"
    
    local basename=$(get_output_basename "$bam_file")
    local success=true
    
    # Run QC tools
    for tool in "${!BAM_QC_OUTPUTS[@]}"; do
        local output_file="${qc_dir}/${basename}${BAM_QC_OUTPUTS[$tool]}"
        
        case "$tool" in
            "FLAGSTAT")
                run_flagstat "$bam_file" "$output_file" || success=false
                ;;
            "QUICKCHECK")
                run_quickcheck "$bam_file" "$output_file" || success=false
                ;;
            "STATS")
                run_stats "$bam_file" "$output_file" || success=false
                ;;
        esac
    done
    
    $success
}

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

#!/bin/bash

source "../config/quality_control_config.sh"

function find_zip_files() {
    local base_dir="$1"
    local pattern="${FILE_PATTERNS[FASTQC_ZIP]}"
    
    log_info "Searching for ZIP files in: $base_dir"
    
    local files=$(find "$base_dir" -type f -name "$pattern")
    
    if [ -z "$files" ]; then
        log_warning "No ZIP files found"
        return 1
    fi
    
    echo "$files"
}

function verify_zip_file() {
    local zip_file="$1"
    
    if ! unzip -t "$zip_file" >/dev/null 2>&1; then
        log_error "Invalid or corrupted ZIP file: $zip_file"
        return 1
    fi
    
    return 0
}

function unzip_in_place() {
    local zip_file="$1"
    local preserve="${2:-${OPERATION_DEFAULTS[PRESERVE_ZIP]}}"
    
    local dir=$(dirname "$zip_file")
    local filename=$(basename "$zip_file")
    
    log_info "Unzipping: $filename"
    
    (
        cd "$dir" || {
            log_error "Failed to change to directory: $dir"
            return 1
        }
        
        if ! unzip -o "$filename"; then
            log_error "Failed to unzip: $filename"
            return 1
        fi
        
        if [ "$preserve" = false ]; then
            rm "$filename" || log_warning "Failed to remove ZIP file: $filename"
        fi
    )
}

function process_zip_files() {
    local files=("$@")
    local total=${#files[@]}
    local success=0
    local failed=0
    
    log_info "Processing $total ZIP files"
    
    for file in "${files[@]}"; do
        if verify_zip_file "$file"; then
            if unzip_in_place "$file"; then
                ((success++))
            else
                ((failed++))
            fi
        else
            ((failed++))
        fi
    done
    
    log_info "Processed files - Success: $success, Failed: $failed"
    return $((failed > 0))
}

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
    fi
    
    if ! samtools index "$output_bam"; then
        log_error "Index creation failed for: $output_bam"
        return 1
    fi
    
    log_info "Completed alignment: $output_bam"
    return 0
}
