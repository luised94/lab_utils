#!/usr/bin/env Rscript
# functions moved

source("../functions/bam_qc_analyzer.R")
source("../functions/mapping_calculator.R")

#' Main BAM QC analysis function
main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    
    if (length(args) == 0) {
        log_error("No directory specified")
        stop("Usage: Rscript analyze_bam_qc.R <directory>")
    }
    
    # Find directory
    directory <- find_data_directory(args[1])
    
    if (is.null(directory)) {
        log_error("Directory not found:", args[1])
        stop("Invalid directory")
    }
    
    # Analyze QC
    results <- analyze_bam_qc(directory)
    
    if (!is.null(results)) {
        print(results)
    }
}

if (!interactive()) {
    main()
}
