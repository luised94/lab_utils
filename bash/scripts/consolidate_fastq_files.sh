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
