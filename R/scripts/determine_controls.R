#!/usr/bin/env Rscript

source("../functions/file_validator.R")
source("../functions/control_manager.R")

#' Main control determination function
determine_sample_controls <- function(directory,
                                    task_id,
                                    reference_pattern = CONFIG$PATTERNS$REFERENCE_GENOME) {
    log_info("Starting control determination")
    
    # Load and validate sample table
    sample_table <- load_sample_table(directory)
    
    # Get matching factors
    factors <- get_factors_to_match(sample_table)
    
    # Setup directories
    bam_dir <- file.path(directory, CONFIG$PATHS$SUBDIRS$ALIGNMENT)
    
    # Get sample BAM
    sample_bam <- find_matching_bam(
        sample_table$sample_ID[task_id],
        bam_dir
    )
    
    if (is.null(sample_bam)) {
        log_error("Sample BAM not found")
        return(NULL)
    }
    
    # Find control
    control_bam <- find_valid_control(
        sample_table[task_id, ],
        sample_table,
        bam_dir,
        factors
    )
    
    if (is.null(control_bam)) {
        log_error("No valid control found")
        return(NULL)
    }
    
    # Return paths
    cat(sample_bam, control_bam, sep = "\n")
}

#' Main entry point
main <- function() {
    if (!interactive()) {
        args <- commandArgs(trailingOnly = TRUE)
        tryCatch({
            arg_list <- validate_input(args)
            determine_sample_controls(
                arg_list$directory_path,
                arg_list$slurm_array_task_id
            )
        }, error = function(e) {
            log_error("Execution failed:", e$message)
            quit(status = 1)
        })
    }
}

if (!interactive()) {
    main()
}
