#!/bin/bash

source "$HOME/lab_utils/bash/functions/slurm_wrapper.sh"
#' Submit SLURM Job Main Function
#' @param array_range Character SLURM array specification
#' @param script_name Character Script name
#' @param dir_name Character Directory name
#' @return Integer 0 if successful
submit_slurm_job_main() {
    if [[ $# -ne 3 ]]; then
        show_usage
        return 1
    fi
    
    # Initialize logging with timestamp for concurrent access
    local timestamp_for_log=$(date +%s)
    local log_file
    log_file=$(initialize_logging "slurm_submit_${timestamp_for_log}")
    
    local array_range="$1"
    local script_name="$2"
    local dir_name="$3"
    local timestamp_for_data=$(date "+%Y-%m-%d-%H-%M-%S")
    
    # Validate inputs
    if ! validate_array_range "$array_range" "$log_file"; then
        return 1
    fi
    
    local script_path
    if ! script_path=$(find_slurm_script "$script_name" "$log_file"); then
        return 1
    fi
    
    local exp_dir
    if ! exp_dir=$(validate_experiment_dir "$dir_name" "$log_file"); then
        return 1
    fi
    
    # Submit job
    log_info "Submitting SLURM job:" "$log_file"
    log_info "  Directory: $exp_dir" "$log_file"
    log_info "  Script: $script_path" "$log_file"
    log_info "  Array: $array_range" "$log_file"
    
    local job_id
    if ! job_id=$(submit_slurm_array_job "$array_range" "$script_name" "$exp_dir" "$job_type"  "$timestamp_for_data" "$log_file"); then
        return 1
    fi
    
    if [[ -z "$job_id" ]]; then
        log_error "Job submission failed" "$log_file"
        return 1
    fi
    
    print_job_info "$job_id"  "$log_file"
    return 0
}


show_usage() {
    cat << EOF
Usage: $(basename "$0") <array_range> <script_name> <directory>

Arguments:
    array_range    SLURM array specification:
                   - Single task: "1"
                   - Multiple tasks: "1,2,5"
                   - Range: "1-10" (automatically adds %16)
    script_name    Script to execute
    directory      Experiment directory name

Examples:
    $(basename "$0") "1-10" "align_fastq.sh" "240304Bel"
    $(basename "$0") "1" "align_fastq.sh" "240304Bel"
    $(basename "$0") "1,2,5" "align_fastq.sh" "240304Bel"

Note: submit_slurm_job automatically applies %16 limit to ranges
EOF
}

# Execute if run as script
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    submit_slurm_job_main "$@"
fi

# bash/scripts/submit_test_job.sh
#!/bin/bash

source "$HOME/lab_utils/bash/functions/slurm_wrapper.sh"

# Calculate array size based on experiment
experiment_dir="241028Bel"
array_range="1-4"  # For testing

# Submit test job
submit_slurm_job_main \
    "$array_range" \
    "test_slurm_settings.sh" \
    "$experiment_dir"

#!/bin/bash
# bash/config/slurm_config.sh

#' SLURM Configuration Settings
#' @description Centralized SLURM configuration for the project

# Load project configuration
source "$HOME/lab_utils/bash/config/project_config.sh"

#' Core SLURM Settings
declare -A SLURM_CONFIG=(
    # Resource Management
    [NODES]="1"
    [TASKS]="1"
    [CPUS_PER_TASK]="4"
    [MEM_PER_CPU]="50G"
    [MAX_ARRAY_SIZE]="16"
    [ARRAY_FORMAT]="1-%d%%16"
    
    # Job Control
    [NICE]="10000"
    [EXCLUDE_NODES]="c[5-22]"
    
    # Notification
    [MAIL_TYPE]="ALL"
    [MAIL_USER]="luised94@mit.edu"
    
    # Paths
    [LOG_DIR]="slurm_logs"
    [MAX_LOG_AGE]="30"
    
)

#' Module Management
declare -A SLURM_MODULES=(
    # Core Modules
    [GNU]="gnu/5.4.0"
    [SAMTOOLS]="samtools/1.10"
    [BOWTIE2]="bowtie2/2.3.5.1"
    [FASTQC]="fastqc/0.11.5"
    
    # Module Loading Order
    [LOAD_ORDER]="GNU BOWTIE2 SAMTOOLS FASTQC"
)

#' Resource Configurations
declare -A SLURM_RESOURCES=(
    # Alignment Settings
    [ALIGN_MEM]="50G"
    [ALIGN_CPUS]="4"
    [ALIGN_TIME]="24:00:00"
    
    # QC Settings
    [QC_MEM]="20G"
    [QC_CPUS]="2"
    [QC_TIME]="4:00:00"
    
    # BigWig Settings
    [BW_MEM]="30G"
    [BW_CPUS]="2"
    [BW_TIME]="8:00:00"
)

#' Reference Genome Settings
declare -A SLURM_GENOMES=(
    [BASE_DIR]="$HOME/data/REFGENS"
    [DEFAULT_GENOME]="SaccharomycescerevisiaeS288C"
    [GENOME_PATTERN]="*_refgenome.fna"
    [INDEX_SUFFIX]="_index"
)

#' Validation Functions
validate_slurm_config() {
    local log_file="$1"
    
    # Verify required directories
    if [[ ! -d "${SLURM_GENOMES[BASE_DIR]}" ]]; then
        log_error "Reference genome directory not found: ${SLURM_GENOMES[BASE_DIR]}" "$log_file"
        return 1
    fi
    
    # Verify module availability
    for module in ${SLURM_MODULES[LOAD_ORDER]}; do
        if ! module avail "${SLURM_MODULES[$module]}" 2>/dev/null; then
            log_error "Required module not available: ${SLURM_MODULES[$module]}" "$log_file"
            return 1
        fi
    done
    
    return 0
}

#' Get SLURM Options
#' @param job_type Character Type of job (ALIGN, QC, BW)
#' @return String SLURM options
get_slurm_options() {
    local job_type="$1"
    
    echo "--nodes=${SLURM_CONFIG[NODES]} \
          --cpus-per-task=${SLURM_RESOURCES[${job_type}_CPUS]} \
          --mem-per-cpu=${SLURM_RESOURCES[${job_type}_MEM]} \
          --time=${SLURM_RESOURCES[${job_type}_TIME]} \
          --exclude=${SLURM_CONFIG[EXCLUDE_NODES]} \
          --nice=${SLURM_CONFIG[NICE]} \
          --mail-type=${SLURM_CONFIG[MAIL_TYPE]} \
          --mail-user=${SLURM_CONFIG[MAIL_USER]}"
}

#!/bin/bash
# bash/functions/slurm_validator.sh

source "$HOME/lab_utils/bash/config/slurm_config.sh"
source "$HOME/lab_utils/bash/functions/logging_utils.sh"

#' Validate SLURM Environment
#' @param job_type Character Type of job (ALIGN, QC, BW)
#' @param log_file Character Log file path
#' @return Integer 0 if valid
validate_slurm_environment() {
    local job_type="$1"
    local log_file="$2"
    
    # Validate resource requirements
    if [[ ! " ALIGN QC BW " =~ " $job_type " ]]; then
        log_error "Invalid job type: $job_type" "$log_file"
        return 1
    fi
    
    # Check SLURM environment
    if [[ -z "${SLURM_JOB_ID:-}" ]]; then
        log_error "Not running in SLURM environment" "$log_file"
        return 1
    fi
    
    # Verify reference directory
    if [[ ! -d "${SLURM_GENOMES[BASE_DIR]}" ]]; then
        log_error "Reference genome directory not found: ${SLURM_GENOMES[BASE_DIR]}" "$log_file"
        return 1
    fi
    
    return 0
}

#' Validate Module Availability
#' @param modules Array Required module names
#' @param log_file Character Log file path
#' @return Integer 0 if all modules available
validate_modules() {
    local -a modules=("$@")
    local log_file="${modules[-1]}" # Last argument is log file
    unset 'modules[-1]'            # Remove log file from array
    
    for module in "${modules[@]}"; do
        if ! module avail "$module" 2>/dev/null; then
            log_error "Required module not available: $module" "$log_file"
            return 1
        fi
    done
    
    return 0
}

#' Format Array Range for MIT SLURM
#' @param range_input String Input range (e.g., "1-10", "1,2,3", "5")
#' @return String Formatted range with %16
format_array_range() {
    local range_input="$1"
    
    # Single number
    if [[ "$range_input" =~ ^[0-9]+$ ]]; then
        echo "$range_input"
        return 0
    fi
    
    # Comma-separated list
    if [[ "$range_input" =~ ^[0-9]+(,[0-9]+)*$ ]]; then
        echo "$range_input"
        return 0
    fi
    
    # Range (add %16 if not present)
    if [[ "$range_input" =~ ^[0-9]+-[0-9]+$ ]]; then
        echo "${range_input}%16"
        return 0
    fi
    
    # Invalid format
    return 1
}
#' Validate SLURM Array Range
#' @param range_input String Input range
#' @param log_file Character Log file path
#' @return Integer 0 if valid
validate_array_range() {
    local range_input="$1"
    local log_file="$2"
    
    # Validate format
    if ! [[ "$range_input" =~ ^([0-9]+|[0-9]+-[0-9]+|[0-9]+(,[0-9]+)*)$ ]]; then
        log_error "Invalid array range format: $range_input" "$log_file"
        return 1
    fi
    
    # Format range
    local formatted_range
    if ! formatted_range=$(format_array_range "$range_input"); then
        log_error "Failed to format array range" "$log_file"
        return 1
    fi
    
    echo "$formatted_range"
    return 0
}


#!/bin/bash
# bash/functions/slurm_wrapper.sh

source "$HOME/lab_utils/bash/config/slurm_config.sh"
source "$HOME/lab_utils/bash/functions/logging_utils.sh"
source "$HOME/lab_utils/bash/functions/slurm_validator.sh"

#' Submit SLURM Array Job
#' @param array_range Character SLURM array specification
#' @param script_name Character Script name
#' @param dir_name Character Directory name
#' @param job_type Character Job type
#' @param log_file Character Log file path
#' @return Integer 0 if successful
submit_slurm_array_job() {
    local array_range="$1"
    local script_name="$2"
    local dir_name="$3"
    local job_type="$4"
    local timestamp_for_data="$5"
    local log_file="$6"
    
    local slurm_opts
    slurm_opts=$(get_slurm_options "$job_type")
    
    local job_id
    job_id=$(sbatch --parsable \
                    --array="$array_range" \
                    $slurm_opts \
                    "$script_name" \
                    "$dir_name"
                    "$timestamp")
    
    if [[ -z "$job_id" ]]; then
        log_error "Job submission failed" "$log_file"
        return 1
    fi
    
    echo "$job_id"
    return 0
}


#' Print Job Information
#' @param job_id Character SLURM job ID
#' @param log_file Character Log file path
print_job_info() {
    local job_id="$1"
    local log_file="$2"
    
    cat << EOF

Job Information:
  Job ID: $job_id
  Log File: $log_file

Monitor Commands:
  squeue -j $job_id
  sacct -j $job_id
  tail -f $log_file
  
Status Check:
  squeue -u $USER
  scontrol show job $job_id
EOF
}

#' Find Script in Project
#' @param script_name Character Script name
#' @param log_file Character Log file path
#' @return String Script path
find_slurm_script() {
    local script_name="$1"
    local log_file="$2"
    
    local script_path="$HOME/lab_utils/bash/scripts/$script_name"
    
    if [[ ! -f "$script_path" ]]; then
        log_error "Script not found: $script_path" "$log_file"
        return 1
    fi
    
    if [[ ! -x "$script_path" ]]; then
        log_error "Script not executable: $script_path" "$log_file"
        return 1
    fi
    
    echo "$script_path"
}

#' Validate Experiment Directory
#' @param dir_name Character Directory name
#' @param log_file Character Log file path
#' @return String Full directory path
validate_experiment_dir() {
    local dir_name="$1"
    local log_file="$2"
    
    local exp_dir="$HOME/data/$dir_name"
    
    if [[ ! -d "$exp_dir" ]]; then
        log_error "Experiment directory not found: $exp_dir" "$log_file"
        return 1
    fi
    
    echo "$exp_dir"
}

#!/bin/bash

source "../config/genome_index_config.sh"

function setup_logging() {
    local log_dir="${GENOME_CONFIG[LOG_DIR]}"
    local timestamp=$(date "+%Y-%m-%d-%H-%M-%S")
    local job_id="${SLURM_ARRAY_JOB_ID:-standalone}"
    
    mkdir -p "$log_dir" || {
        echo "Failed to create log directory: $log_dir"
        return 1
    }
    
    echo "${log_dir}/indexing_${job_id}_${timestamp}"
}

function validate_slurm_env() {
    if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
        log_error "This script must be run as a SLURM array job"
        return 1
    fi
    
    if [ -z "${SLURM_ARRAY_JOB_ID:-}" ]; then
        log_error "SLURM_ARRAY_JOB_ID not set"
        return 1
    fi
    
    return 0
}

function load_required_modules() {
    log_info "Loading required modules"
    
    module purge || log_warning "Failed to purge modules"
    
    for module in "${MODULES[@]}"; do
        if ! module load "$module"; then
            log_error "Failed to load module: $module"
            return 1
        fi
        log_info "Loaded module: $module"
    done
    
    return 0
}

#!/bin/bash

source "../config/slurm_config.sh"
# Assuming logging utilities are already sourced

function find_slurm_files() {
    local search_dir="${1:-.}"
    local pattern="${2:-$SLURM_OUTPUT_PATTERN}"
    local max_depth="${3:-$MAX_DEPTH_SEARCH}"

    log_info "Searching for SLURM files in: $search_dir"
    
    if [[ ! -d "$search_dir" ]]; then
        log_error "Directory does not exist: $search_dir"
        return 1
    fi

    find "$search_dir" \
        -maxdepth "$max_depth" \
        -type f \
        -name "$pattern" \
        -printf "%T@ %p\n" | \
        sort -nr | \
        cut -d' ' -f2-
}

function organize_slurm_files() {
    local source_dir="${1:-.}"
    local target_dir="${2:-$DEFAULT_SLURM_LOG_DIR}"
    
    log_info "Organizing SLURM files from $source_dir to $target_dir"
    
    # Create target directory if it doesn't exist
    mkdir -p "$target_dir"
    
    # Move files to organized structure
    while IFS= read -r file; do
        local job_id=$(basename "$file" | sed 's/slurm-\([0-9]*\).out/\1/')
        local date_str=$(stat -c %y "$file" | cut -d' ' -f1)
        local year_month=$(date -d "$date_str" +%Y/%m)
        local target_subdir="$target_dir/$year_month"
        
        mkdir -p "$target_subdir"
        mv "$file" "$target_subdir/"
        log_info "Moved $file to $target_subdir/"
    done < <(find_slurm_files "$source_dir")
}

function cleanup_old_logs() {
    local log_dir="${1:-$DEFAULT_SLURM_LOG_DIR}"
    local max_age="${2:-$MAX_LOG_AGE_DAYS}"
    
    log_info "Cleaning up SLURM logs older than $max_age days"
    
    find "$log_dir" \
        -type f \
        -name "$SLURM_OUTPUT_PATTERN" \
        -mtime +"$max_age" \
        -exec rm {} \;
}

function check_large_files() {
    local search_dir="${1:-.}"
    local max_size="${2:-$MAX_FILE_SIZE_MB}"
    
    log_info "Checking for large SLURM output files"
    
    find "$search_dir" \
        -type f \
        -name "$SLURM_OUTPUT_PATTERN" \
        -size +"${max_size}M" \
        -exec ls -lh {} \; | \
    while read -r line; do
        log_warning "Large SLURM output file: $line"
    done
}

#!/bin/bash
# functions moved

source "../functions/slurm_output_manager.sh"

function show_usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]
Manage SLURM output files

Options:
    -d, --directory DIR    Search in specific directory
    -o, --organize        Organize files into dated structure
    -c, --cleanup        Remove old log files
    -l, --large          Check for large files
    -h, --help           Show this help message
EOF
}

function main() {
    local search_dir="."
    local do_organize=false
    local do_cleanup=false
    local check_large=false

    while [[ $# -gt 0 ]]; do
        case $1 in
            -d|--directory) search_dir="$2"; shift 2 ;;
            -o|--organize) do_organize=true; shift ;;
            -c|--cleanup) do_cleanup=true; shift ;;
            -l|--large) check_large=true; shift ;;
            -h|--help) show_usage; exit 0 ;;
            *) log_error "Unknown option: $1"; show_usage; exit 1 ;;
        esac
    done

    if $do_organize; then
        organize_slurm_files "$search_dir"
    fi

    if $do_cleanup; then
        cleanup_old_logs
    fi

    if $check_large; then
        check_large_files "$search_dir"
    fi

    if ! $do_organize && ! $do_cleanup && ! $check_large; then
        find_slurm_files "$search_dir"
    fi
}

main "$@"

#!/bin/bash

source "../config/slurm_config.sh"

function validate_array_range() {
    local range="$1"
    
    log_info "Validating array range: $range"
    
    if [[ ! "$range" =~ ^[0-9]+(|-[0-9]+(%[0-9]+)?|,[0-9]+)*$ ]]; then
        log_error "Invalid array range format: $range"
        return 1
    fi
    
    if [[ "$range" =~ %([0-9]+) ]]; then
        local concurrent="${BASH_REMATCH[1]}"
        if ((concurrent > ${SLURM_WRAPPER[MAX_ARRAY_SIZE]})); then
            log_warning "Concurrent jobs ($concurrent) exceeds recommended maximum (${SLURM_WRAPPER[MAX_ARRAY_SIZE]})"
        }
        fi
    
    return 0
}

function find_script() {
    local script_name="$1"
    local base_dir="${SLURM_WRAPPER[SCRIPT_BASE_DIR]}"
    
    log_info "Searching for script: $script_name"
    
    local script_path=$(find "$base_dir" -type f -name "$script_name")
    
    if [ -z "$script_path" ]; then
        log_error "Script not found: $script_name"
        return 1
    fi
    
    if [ ! -x "$script_path" ]; then
        log_error "Script not executable: $script_path"
        return 1
    fi
    
    echo "$script_path"
}

function validate_experiment_dir() {
    local dir_name="$1"
    local base_dir="${SLURM_WRAPPER[DATA_BASE_DIR]}"
    
    log_info "Validating experiment directory: $dir_name"
    
    local full_path=$(find -H "$base_dir" -maxdepth 1 -type d -name "$dir_name")
    
    if [ -z "$full_path" ]; then
        log_error "Directory not found: $dir_name"
        return 1
    fi
    
    echo "$full_path"
}

#!/bin/bash
# bash/scripts/test_slurm_settings.sh

source "$HOME/lab_utils/bash/functions/logging_utils.sh"

#' Test SLURM Settings
#' @param experiment_dir Character Experiment directory
#' @return Integer 0 if successful
test_slurm_settings_main() {
    if [[ $# -ne 1 ]]; then
        echo "Usage: $0 <experiment_dir>"
        return 1
    }
    
    local experiment_dir="$1"
    local log_file
    log_file=$(initialize_logging "test_slurm")
    
    # Print SLURM Environment
    cat << EOF
=== SLURM Job Information ===
Job ID: ${SLURM_JOB_ID:-Not Set}
Array Job ID: ${SLURM_ARRAY_JOB_ID:-Not Set}
Array Task ID: ${SLURM_ARRAY_TASK_ID:-Not Set}
Job Name: ${SLURM_JOB_NAME:-Not Set}

=== SLURM Resource Allocation ===
Nodes: ${SLURM_JOB_NODELIST:-Not Set}
Node Count: ${SLURM_NNODES:-Not Set}
CPUs per Task: ${SLURM_CPUS_PER_TASK:-Not Set}
Tasks per Node: ${SLURM_NTASKS_PER_NODE:-Not Set}
Memory per Node: ${SLURM_MEM_PER_NODE:-Not Set}

=== Directory Information ===
Working Directory: ${SLURM_SUBMIT_DIR:-Not Set}
Experiment Directory: $experiment_dir

=== Module Status ===
$(module list 2>&1)

=== System Information ===
Hostname: $(hostname)
Current Directory: $(pwd)
User: $USER
Date: $(date)
EOF

    return 0
}

# Execute if run as script
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    test_slurm_settings_main "$@"
fi

#!bin/bash

find . -maxdepth 1 -type f -name "slurm-*.out" -exec echo {} +
