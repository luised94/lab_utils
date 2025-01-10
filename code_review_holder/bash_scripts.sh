
#!/bin/bash
# functions moved

set -euo pipefail

source "../functions/archive_handler.sh"

function show_usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS] <directory>
Unzips FASTQC result files in specified directory

Options:
    -f, --force       Skip confirmation
    -k, --keep        Keep ZIP files after extraction
    -b, --batch SIZE  Process in batches of SIZE files
    -h, --help        Show this help message

Example:
    $(basename "$0") experiment_20240101
EOF
}

function confirm_operation() {
    local files="$1"
    local timeout="${OPERATION_DEFAULTS[CONFIRM_TIMEOUT]}"
    
    echo "Files to process:"
    echo "$files"
    echo
    
    read -t "$timeout" -p "Proceed with unzipping these files? (y/n): " -r || {
        echo
        log_error "Confirmation timed out"
        return 1
    }
    
    [[ $REPLY =~ ^[Yy]$ ]]
}

function validate_directory() {
    local dir_name="$1"
    local base_dir="${QC_PATHS[BASE_DIR]}"
    
    log_info "Validating directory: $dir_name"
    
    local full_path=$(find -H "$base_dir" \
                      -maxdepth "${QC_PATHS[MAX_DEPTH]}" \
                      -type d \
                      -name "$dir_name")
    
    if [ -z "$full_path" ]; then
        log_error "Directory not found: $dir_name"
        return 1
    }
    
    echo "$full_path"
}

function main() {
    local force=false
    local preserve=${OPERATION_DEFAULTS[PRESERVE_ZIP]}
    local batch_size=${OPERATION_DEFAULTS[UNZIP_BATCH_SIZE]}
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -f|--force) force=true; shift ;;
            -k|--keep) preserve=true; shift ;;
            -b|--batch) batch_size="$2"; shift 2 ;;
            -h|--help) show_usage; exit 0 ;;
            *) break ;;
        esac
    done
    
    if [ $# -ne 1 ]; then
        log_error "Directory name required"
        show_usage
        exit 1
    fi
    
    local exp_dir=$(validate_directory "$1") || exit 1
    local qc_dir="$exp_dir/${QC_PATHS[QC_SUBDIR]}"
    
    mapfile -t zip_files < <(find_zip_files "$qc_dir") || {
        log_error "No ZIP files found in: $qc_dir"
        exit 1
    }
    
    if [ "$force" = false ]; then
        confirm_operation "${zip_files[*]}" || exit 1
    fi
    
    process_zip_files "${zip_files[@]}"
}

main "$@"

#!/bin/bash
# bash/scripts/transfer_bmc_experiment_to_luria.sh
# Requires password twice, for the transfer and verification.
# Run after setup_bmc_experiment.R

# Source dependencies
source "$HOME/lab_utils/bash/functions/logging_utils.sh"

#' Check Network Connectivity
#' @param remote_host Character Remote host address
#' @return Integer 0 if successful, 1 otherwise
check_connection() {
    local remote_host="$1"
    
    if ! ping -c 1 "$remote_host" &> /dev/null; then
        log_error "Cannot connect to $remote_host"
        log_info "Please ensure:"
        log_info "1. VPN is connected"
        log_info "2. You can connect via: ssh -A -Y ${PROJECT_CONFIG[REMOTE_USER]}@$remote_host"
        return 1
    fi
    return 0
}

#' Validate Local Directory Structure
#' @param dir Character Directory to validate
#' @return Integer 0 if valid, 1 otherwise
validate_directory() {
    local dir="$1"
    local log_file="$2"
    
    if [ ! -d "$dir" ]; then
        log_error "Directory not found: $dir" "${log_file}"
        log_info "Provide full path." "${log_file}"
        return 1
    fi
    
    # Check required subdirectories
    for subdir in ${PROJECT_CONFIG[REQUIRED_DIRS]}; do
        if [ ! -d "$dir/$subdir" ]; then
            log_error "Required subdirectory missing: $subdir"
            return 1
        fi
    done
    return 0
}

#' Transfer Data to Remote Host
#' @param source_dir Character Source directory
#' @param log_file Character Log file path
#' @return Integer 0 if successful, 1 otherwise
transfer_data() {
    local source_dir="$1"
    local log_file="$2"
    local experiment_id=$(basename "$source_dir")
    
    log_info "Starting transfer of $experiment_id" "$log_file"
    
    rsync -avzP --stats \
        "$source_dir/" \
        "${PROJECT_CONFIG[REMOTE_USER]}@${PROJECT_CONFIG[REMOTE_HOST]}:${PROJECT_CONFIG[REMOTE_PATH]}/$experiment_id/" \
        2>&1 | tee -a "$log_file"
    
    if [ ${PIPESTATUS[0]} -eq 0 ]; then
        log_info "Transfer completed successfully" "$log_file"
        return 0
    else
        log_error "Transfer failed" "$log_file"
        return 1
    fi
}

#' Verify Transfer Completion
#' @param source_dir Character Source directory
#' @param log_file Character Log file path
#' @return Integer 0 if verified, 1 otherwise
verify_transfer() {
    local source_dir="$1"
    local log_file="$2"
    local experiment_id=$(basename "$source_dir")
    
    log_info "Verifying transfer" "$log_file"
    
    local local_count=$(find "$source_dir" -type f | wc -l)
    local remote_count=$(ssh "${PROJECT_CONFIG[REMOTE_USER]}@${PROJECT_CONFIG[REMOTE_HOST]}" \
        "find ${PROJECT_CONFIG[REMOTE_PATH]}/$experiment_id -type f | wc -l")
    
    if [ "$local_count" -eq "$remote_count" ]; then
        log_info "Verification successful: $local_count files transferred" "$log_file"
        return 0
    else
        log_error "Verification failed: Local=$local_count Remote=$remote_count" "$log_file"
        return 1
    fi
}

#' Main Function
#' @param args Array Script arguments
#' @return None
transfer_bmc_experiment_to_luria_main() {
    # Initialize logging
    local log_file
    log_file="$(initialize_logging "transfer_bmc_experiment")"
    
    if [ $# -ne 1 ]; then
        log_error "Usage: $0 <experiment_directory>" "$log_file"
        exit 1
    fi
    
    local source_dir="$1"
    
    # Run checks
    check_connection "${PROJECT_CONFIG[REMOTE_HOST]}" || exit 1
    validate_directory "$source_dir" "$log_file" || exit 1
    
    # Transfer and verify
    transfer_data "$source_dir" "$log_file" || exit 1
    verify_transfer "$source_dir" "$log_file" || exit 1
    
    log_info "Next steps:" "$log_file"
    log_info "1. Login to cluster: ssh -A -Y ${PROJECT_CONFIG[REMOTE_USER]}@${PROJECT_CONFIG[REMOTE_HOST]}" "$log_file"
    log_info "2. Run: bash ~/lab_utils/bash/scripts/download_bmc_data.sh" "$log_file"
}

# Execute if run as script
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    transfer_bmc_experiment_to_luria_main "$@"
fi

#!/bin/bash
# functions moved

#SBATCH parameters from SLURM_CONFIG

set -euo pipefail

source "../functions/quality_control.sh"

function validate_input() {
    local dir_name="$1"
    
    if [ -z "$dir_name" ]; then
        log_error "Directory name not provided"
        return 1
    }
    
    local full_path="$HOME/data/$dir_name"
    if [ ! -d "$full_path" ]; then
        log_error "Directory not found: $full_path"
        return 1
    }
    
    echo "$full_path"
}

function main() {
    if [ $# -lt 1 ]; then
        log_error "Usage: $0 <directory_name>"
        exit 1
    }
    
    local base_dir=$(validate_input "$1") || exit 1
    local qc_dir=$(setup_qc_directory "$base_dir") || exit 1
    
    load_required_modules || exit 1
    
    mapfile -t raw_files < <(find_fastq_files "$base_dir")
    mapfile -t processed_files < <(find_processed_fastq "$base_dir")
    
    local task_index=$((SLURM_ARRAY_TASK_ID - 1))
    
    if [ -n "${raw_files[$task_index]:-}" ]; then
        run_fastqc "${raw_files[$task_index]}" "$qc_dir" || exit 1
    fi
    
    if [ -n "${processed_files[$task_index]:-}" ]; then
        run_fastqc "${processed_files[$task_index]}" "$qc_dir" || exit 1
    fi
    
    log_info "Quality control completed successfully"
}

main "$@"

#!/bin/bash
# functions moved

#SBATCH parameters from SLURM_CONFIG

set -euo pipefail

source "../functions/fastq_processor.sh"

function main() {
    if [ $# -ne 2 ]; then
        log_error "Usage: $0 <directory> <time_id>"
        exit 1
    }
    
    local exp_dir="$HOME/data/$1"
    local time_id="$2"
    
    log_info "Starting FASTQ processing"
    log_info "Job ID: ${SLURM_JOB_ID}"
    log_info "Array Task: ${SLURM_ARRAY_TASK_ID}"
    
    setup_output_directories "$exp_dir" || exit 1
    
    # Load required modules
    module purge
    for module in "${REQUIRED_MODULES[@]}"; do
        if ! module load "$module"; then
            log_error "Failed to load module: $module"
            exit 1
        fi
    done
    
    # Get input files
    mapfile -t input_files < <(find_input_files "$exp_dir")
    
    if [ ${#input_files[@]} -eq 0 ]; then
        log_error "No input files found"
        exit 1
    }
    
    local task_index=$((SLURM_ARRAY_TASK_ID - 1))
    if [ $task_index -ge ${#input_files[@]} ]; then
        log_error "Task ID exceeds number of files"
        exit 1
    }
    
    process_fastq_file "${input_files[$task_index]}" "$exp_dir" || exit 1
    
    log_info "Processing completed successfully"
}

main "$@"
#!/bin/bash
# functions moved

#SBATCH parameters from SLURM_CONFIG

set -euo pipefail

source "../functions/bam_processor.sh"

function main() {
    if [ $# -ne 2 ]; then
        log_error "Usage: $0 <directory> <time_id>"
        exit 1
    }
    
    local exp_dir="$HOME/data/$1"
    local time_id="$2"
    
    log_info "Starting BAM quality control"
    log_info "Job ID: ${SLURM_JOB_ID}"
    log_info "Array Task: ${SLURM_ARRAY_TASK_ID}"
    
    setup_qc_directories "$exp_dir" || exit 1
    
    # Load required modules
    module purge
    for module in "${REQUIRED_MODULES[@]}"; do
        if ! module load "$module"; then
            log_error "Failed to load module: $module"
            exit 1
        fi
    done
    
    # Get BAM files
    mapfile -t bam_files < <(find_bam_files "$exp_dir")
    
    if [ ${#bam_files[@]} -eq 0 ]; then
        log_error "No BAM files found"
        exit 1
    }
    
    local task_index=$((SLURM_ARRAY_TASK_ID - 1))
    if [ $task_index -ge ${#bam_files[@]} ]; then
        log_error "Task ID exceeds number of files"
        exit 1
    }
    
    process_bam_file "${bam_files[$task_index]}" "${exp_dir}/${QC_DIRS[OUTPUT]}" || exit 1
    
    log_info "Quality control completed successfully"
}

main "$@"

#!/bin/bash
# functions moved

#SBATCH parameters from SLURM_CONFIG

set -euo pipefail

source "../functions/r_integration.sh"
source "../functions/bam_comparer.sh"

function main() {
    if [ $# -ne 2 ]; then
        log_error "Usage: $0 <directory> <time_id>"
        exit 1
    }
    
    local exp_dir="$HOME/data/$1"
    local time_id="$2"
    
    log_info "Starting BAM comparison"
    log_info "Job ID: ${SLURM_JOB_ID}"
    log_info "Array Task: ${SLURM_ARRAY_TASK_ID}"
    
    # Load required modules
    module purge
    for module in "${REQUIRED_MODULES[@]}"; do
        if ! module load "$module"; then
            log_error "Failed to load module: $module"
            exit 1
        fi
    done
    
    # Get BAM pairs from R
    mapfile -t bam_pairs < <(get_bam_pairs "$(basename "$exp_dir")" "$SLURM_ARRAY_TASK_ID")
    
    validate_bam_pairs "${bam_pairs[@]}" || exit 1
    
    local output_file=$(generate_output_name "${bam_pairs[0]}" \
                                           "${bam_pairs[1]}" \
                                           "$time_id" \
                                           "${exp_dir}/${OUTPUT_DIRS[BIGWIG]}")
    
    local threads=$((SLURM_CPUS_PER_TASK / 2))
    
    run_comparison "${bam_pairs[0]}" "${bam_pairs[1]}" "$output_file" "$threads" || exit 1
    
    log_info "Comparison completed successfully"
}

main "$@"


# functions moved

#SBATCH parameters from SLURM_CONFIG

set -euo pipefail

source "../functions/alignment_handler.sh"

function setup_directories() {
    local base_dir="$1"
    
    local dirs=(
        "logs"
        "alignment"
    )
    
    for dir in "${dirs[@]}"; do
        local full_path="${base_dir}/${dir}"
        mkdir -p "$full_path" || {
            log_error "Failed to create directory: $full_path"
            return 1
        }
    done
}

function main() {
    if [ $# -ne 2 ]; then
        log_error "Usage: $0 <experiment_dir> <time_id>"
        exit 1
    }
    
    local exp_dir="$HOME/data/$1"
    local time_id="$2"
    
    setup_directories "$exp_dir" || exit 1
    
    log_info "Starting alignment process"
    log_info "Job ID: ${SLURM_JOB_ID}"
    log_info "Array Task: ${SLURM_ARRAY_TASK_ID}"
    
    # Load required modules
    module purge
    for module in "${MODULE_REQUIREMENTS[@]}"; do
        if ! module load "$module"; then
            log_error "Failed to load module: $module"
            exit 1
        fi
    done
    
    # Get file lists
    mapfile -t fastq_files < <(find "$exp_dir" -type f -name "${FILE_PATTERNS[FASTQ]}")
    mapfile -t genome_files < <(find "$REFGENOME_DIR" -type f -name "${FILE_PATTERNS[GENOME]}")
    
    if [ ${#fastq_files[@]} -eq 0 ] || [ ${#genome_files[@]} -eq 0 ]; then
        log_error "No input files found"
        exit 1
    }
    
    log_info "Found ${#fastq_files[@]} FASTQ files and ${#genome_files[@]} genomes"
    
    # Calculate indices
    local indices=$(calculate_indices "$SLURM_ARRAY_TASK_ID" "${#fastq_files[@]}")
    local genome_index=${indices%:*}
    local fastq_index=${indices#*:}
    
    perform_alignment "${genome_files[$genome_index]}" \
                     "${fastq_files[$fastq_index]}" \
                     "${exp_dir}/alignment" \
                     "$SLURM_CPUS_PER_TASK" || exit 1
    
    log_info "Alignment task completed successfully"
}

main "$@"

#!/bin/bash
# functions moved

set -euo pipefail

source "../functions/genome_organizer.sh"

function show_usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]
Reorganizes reference genome directories

Options:
    -d, --dir DIR     Process specific directory
    -a, --all         Process all genome directories
    -s, --standardize Standardize chromosome names
    -h, --help        Show this help message
EOF
}

function process_genome_directory() {
    local dir="$1"
    
    log_info "Processing directory: $dir"
    
    local assembly_report="${dir}/${GENOME_PATHS[NCBI_DATA]}/${GENOME_PATHS[ASSEMBLY_REPORT]}"
    
    local organism_name=$(extract_organism_name "$assembly_report") || return 1
    
    if ! reorganize_genome_files "$dir" "$organism_name"; then
        log_error "Failed to reorganize: $dir"
        return 1
    }
    
    if [ -d "$dir" ] && [ "$dir" != "$organism_name" ]; then
        mv "$dir" "$organism_name" || {
            log_error "Failed to rename directory to: $organism_name"
            return 1
        }
    }
    
    log_info "Successfully processed: $organism_name"
    return 0
}

function main() {
    local process_all=false
    local standardize=false
    local target_dir=""
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -d|--dir) target_dir="$2"; shift 2 ;;
            -a|--all) process_all=true; shift ;;
            -s|--standardize) standardize=true; shift ;;
            -h|--help) show_usage; exit 0 ;;
            *) log_error "Unknown option: $1"; show_usage; exit 1 ;;
        esac
    done
    
    if $process_all; then
        mapfile -t dirs < <(find . -maxdepth 1 -type d -name "GC[AF]_*")
    elif [ -n "$target_dir" ]; then
        dirs=("$target_dir")
    else
        log_error "No directory specified"
        show_usage
        exit 1
    fi
    
    for dir in "${dirs[@]}"; do
        process_genome_directory "$dir" || continue
    done
    
    log_info "Genome reorganization completed"
}

main "$@"

#!/bin/bash
# functions moved

set -euo pipefail

source "../functions/fasta_processor.sh"

function show_usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]
Reformats S288C genome headers to UCSC format

Options:
    -m, --method METHOD  Use specific method (awk|while)
    -n, --no-backup     Skip backup creation
    -r, --restore       Restore from backup
    -h, --help          Show this help message
EOF
}

function find_s288c_genome() {
    local base_dir="$1"
    local pattern="${GENOME_FILES[S288C_PATTERN]}"
    
    log_info "Searching for S288C genome"
    
    local genome_path=$(find "$base_dir" -type f -name "$pattern")
    
    if [ -z "$genome_path" ]; then
        log_error "S288C genome not found"
        return 1
    fi
    
    echo "$genome_path"
}

function restore_from_backup() {
    local genome_path="$1"
    local backup_file="${genome_path%_refgenome.fna}${GENOME_FILES[BACKUP_SUFFIX]}"
    
    if [ ! -f "$backup_file" ]; then
        log_error "Backup file not found: $backup_file"
        return 1
    }
    
    log_info "Restoring from backup"
    cp "$backup_file" "$genome_path"
}

function main() {
    local method="awk"
    local do_backup=${FILE_OPERATIONS[BACKUP]}
    local restore=false
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -m|--method) method="$2"; shift 2 ;;
            -n|--no-backup) do_backup=false; shift ;;
            -r|--restore) restore=true; shift ;;
            -h|--help) show_usage; exit 0 ;;
            *) log_error "Unknown option: $1"; show_usage; exit 1 ;;
        esac
    done
    
    local genome_path=$(find_s288c_genome "$REFGENOME_DIR") || exit 1
    
    if [ "$restore" = true ]; then
        restore_from_backup "$genome_path"
        exit 0
    fi
    
    validate_fasta "$genome_path" || exit 1
    
    if [ "$do_backup" = true ]; then
        create_backup "$genome_path" || exit 1
    fi
    
    local backup_file="${genome_path%_refgenome.fna}${GENOME_FILES[BACKUP_SUFFIX]}"
    
    if [ "$method" = "awk" ]; then
        reformat_headers_awk "$backup_file" "$genome_path"
    else
        reformat_headers_while "$backup_file" "$genome_path"
    fi
    
    if [ "${FILE_OPERATIONS[VERIFY]}" = true ]; then
        verify_conversion "$backup_file" "$genome_path" || exit 1
    fi
    
    log_info "Header reformatting completed successfully"
}

main "$@"

#!/usr/bin/env bash
# functions moved

set -o errexit
set -o nounset
set -o pipefail

source "../functions/ngs_file_manager.sh"

function show_usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]
Manage NGS data files

Options:
    -t, --target DIR     Target directory for file movement
    -s, --source DIR     Source directory to search
    -d, --depth NUM      Maximum search depth (default: ${DEFAULTS[MAX_DEPTH]})
    -b, --batch NUM      Batch size for processing (default: ${DEFAULTS[BATCH_SIZE]})
    -a, --analyze        Only analyze file distribution
    -h, --help          Show this help message
EOF
}

function main() {
    local target_dir=""
    local search_dir="."
    local do_analyze=false
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -t|--target) target_dir="$2"; shift 2 ;;
            -s|--source) search_dir="$2"; shift 2 ;;
            -d|--depth) DEFAULTS[MAX_DEPTH]="$2"; shift 2 ;;
            -b|--batch) DEFAULTS[BATCH_SIZE]="$2"; shift 2 ;;
            -a|--analyze) do_analyze=true; shift ;;
            -h|--help) show_usage; exit 0 ;;
            *) log_error "Unknown option: $1"; show_usage; exit 1 ;;
        esac
    done
    
    if $do_analyze; then
        analyze_file_distribution "$search_dir"
    else
        if [[ -z "$target_dir" ]]; then
            log_error "Target directory must be specified"
            show_usage
            exit 1
        fi
        move_ngs_files "$target_dir" "$search_dir"
    fi
}

main "$@"

#!/bin/bash
# functions moved

#SBATCH parameters from SLURM_CONFIG

set -euo pipefail

source "../functions/coverage_processor.sh"

function main() {
    if [ $# -ne 2 ]; then
        log_error "Usage: $0 <directory> <time_id>"
        exit 1
    }
    
    local exp_dir="$HOME/data/$1"
    local time_id="$2"
    
    log_info "Starting coverage generation"
    log_info "Job ID: ${SLURM_JOB_ID}"
    log_info "Array Task: ${SLURM_ARRAY_TASK_ID}"
    
    setup_coverage_directories "$exp_dir" || exit 1
    
    # Load required modules
    module purge
    for module in "${REQUIRED_MODULES[@]}"; do
        if ! module load "$module"; then
            log_error "Failed to load module: $module"
            exit 1
        fi
    done
    
    # Get BAM files
    mapfile -t bam_files < <(find_bam_files "$exp_dir")
    
    if [ ${#bam_files[@]} -eq 0 ]; then
        log_error "No BAM files found"
        exit 1
    }
    
    local task_index=$((SLURM_ARRAY_TASK_ID - 1))
    if [ $task_index -ge ${#bam_files[@]} ]; then
        log_error "Task ID exceeds number of files"
        exit 1
    }
    
    local output_file=$(generate_output_name "${bam_files[$task_index]}" \
                                           "$time_id" \
                                           "${exp_dir}/${OUTPUT_DIRS[BIGWIG]}")
    
    process_bam_coverage "${bam_files[$task_index]}" "$output_file" || exit 1
    
    log_info "Coverage generation completed successfully"
}

main "$@"

#!/bin/bash
# functions moved

set -euo pipefail

source "../functions/sra_downloader.sh"

function show_usage() {
    cat << EOF
Usage: $(basename "$0") <directory>
Downloads Eaton 2010 paper data for control and comparison

Arguments:
    directory    Target directory name (relative to ${SRA_CONFIG[DATA_DIR]})

Study Information:
    BioProject: ${EATON_2010[BIOPROJECT]}
    Description: ${EATON_2010[DESCRIPTION]}
EOF
}

function main() {
    if [ $# -ne 1 ]; then
        show_usage
        exit 1
    }
    
    local output_dir=$(validate_input "$1") || exit 1
    local combined_output="${output_dir}/${EATON_2010[OUTPUT_FILE]}.gz"
    local downloaded_files=()
    
    for file_name in "${!SAMPLES[@]}"; do
        local accession="${SAMPLES[$file_name]}"
        local url=$(construct_download_url "$accession")
        local output_file="${output_dir}/$file_name"
        
        verify_url "$url" || continue
        download_file "$url" "$output_file" || continue
        
        downloaded_files+=("$file_name")
    done
    
    if [ ${#downloaded_files[@]} -gt 0 ]; then
        concatenate_files "$output_dir" "$combined_output" "${downloaded_files[@]}" || exit 1
        decompress_file "$combined_output" || exit 1
    else
        log_error "No files were downloaded successfully"
        exit 1
    fi
    
    log_info "Processing complete"
}

main "$@"

#!/bin/bash
# bash/scripts/fastq/download_bmc_data.sh

#' Download BMC Data Script
#' @description Download and process BMC FASTQ data

echo "Initializing BMC data download"
echo "Script location: ${BASH_SOURCE[0]}"
echo "Working directory: $(pwd)"

# Find repository root using git
repo_root=$(git rev-parse --show-toplevel 2>/dev/null) || {
    echo "? Not in a git repository"
    exit 1
}

# Initialize environment
source "$repo_root/bash/core/initialize_lab_environment.sh" || {
    echo "? Failed to initialize environment"
    exit 1
}


# Load required modules
echo "Loading required modules"
for module in "fastq/bmc_handler" "fastq/fastq_processor"; do
    echo "Loading: $module"
    if ! load_lab_module "$module"; then
        log_error "Failed to load module: $module"
        exit 1
    fi
    echo "Loaded successfully"
done
#' Download BMC Data Main Function
#' @param experiment_id Character Experiment identifier
#' @return Integer 0 if successful
download_bmc_data_main() {
    local experiment_id="$1"
    local log_file
    
    # Initialize logging
    log_file=$(initialize_logging "download_bmc_data") || {
        echo "? Failed to initialize logging"
        return 1
    }
    
    log_info "Starting BMC data download" "$log_file"
    log_debug "Environment verification" "$log_file"
    log_debug "LAB_UTILS_ROOT: $LAB_UTILS_ROOT" "$log_file"
    log_debug "Experiment ID: $experiment_id" "$log_file"
    log_debug "Log file: $log_file" "$log_file"
    
    # Verify host
    if ! verify_host; then
        log_error "Host verification failed" "$log_file"
        return 1
    fi
    
    # Validate paths
    local paths
    log_debug "Validating BMC paths" "$log_file"
    if ! paths=$(validate_bmc_paths "$experiment_id" "$log_file"); then
        return 1
    fi
    
    local bmc_path=${paths%:*}
    local local_path=${paths#*:}

    log_debug " BMC path: $bmc_path" "$log_file"
    log_debug " Local path: $local_path" "$log_file"
    
    # Download data
    if ! download_from_bmc "$paths" "$log_file"; then
        return 1
    fi

    # Organize files
    if ! move_fastq_files_to_current_directory "${local_path}" "$log_file"; then
        return 1
    fi

    # Cleanup
    if ! clean_experiment_directory "${local_path}" "$log_file"; then
        return 1
    fi

    log_info "Download process completed successfully" "$log_file"
    return 0
}

# Show usage information
show_usage() {
    cat << EOF
Usage: $(basename "$0") <experiment_id>
Arguments:
    experiment_id    Experiment identifier (format: YYMMDD'Bel')
Example:
    $(basename "$0") 241010Bel
EOF
}

# Execute if run as script
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    download_bmc_data_main "$@"
fi

#!/bin/bash
# bash/scripts/consolidate_fastq_files.sh

source "$HOME/lab_utils/bash/config/project_config.sh"
source "$HOME/lab_utils/bash/functions/fastq_consolidator_processor.sh"

#' Consolidate FASTQ Files Main Function
#' @param experiment_id Character Experiment identifier
#' @return Integer 0 if successful
consolidate_fastq_files_main() {
    local experiment_id="$1"
    if [[ $# -ne 1 ]]; then
        show_usage
        return 1
    fi
    
    # Initialize logging
    local log_file
    log_file=$(initialize_logging "consolidate_fastq")
    
    local experiment_dir
    
    # Validate input
    if ! experiment_dir=$(validate_fastq_input "$experiment_id" "$log_file"); then
        return 1
    fi
    
    # Setup directories
    local output_dir
    if ! output_dir=$(setup_fastq_directories "$experiment_dir" "$log_file"); then
        return 1
    fi
    
    # Process files
    if ! consolidate_fastq_files_by_id "$experiment_dir" "$output_dir" "$log_file"; then
        return 1
    fi
    
    log_info "FASTQ consolidation completed successfully" "$log_file"
    return 0
}

show_usage() {
    cat << EOF
Usage: $(basename "$0") <experiment_id>

Arguments:
    experiment_id    Experiment identifier (format: YYMMDD'Bel')

Example:
    $(basename "$0") 241028Bel
EOF
}

# Execute if run as script
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    consolidate_fastq_files_main "$@"
fi
