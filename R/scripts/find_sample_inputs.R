#!/usr/bin/env Rscript

source("../functions/sample_matcher.R")
source("../functions/bam_finder.R")

#' Main sample input determination function
determine_input_for_all_samples <- function(sample_table,
                                          directory_path,
                                          reference_pattern = CONFIG$FILES$PATTERNS$REFERENCE) {
    log_info("Starting sample input determination")
    
    # Validate inputs
    validate_inputs(sample_table, directory_path)
    
    # Find matching samples
    results <- find_matching_samples(
        sample_table,
        directory_path
    )
    
    # Output results
    output_results(results)
    
    log_info("Sample input determination completed")
}

validate_inputs <- function(sample_table, directory_path) {
    if (!all(CONFIG$VALIDATION$REQUIRED_COLUMNS %in% colnames(sample_table))) {
        log_error("Missing required columns in sample table")
        stop("Invalid sample table")
    }
    
    if (!dir.exists(directory_path)) {
        log_error("Directory not found:", directory_path)
        stop("Invalid directory")
    }
}

output_results <- function(results) {
    for (result in results) {
        if (!is.null(result)) {
            cat(result$sample, result$control, sep = "\n")
        }
    }
}
