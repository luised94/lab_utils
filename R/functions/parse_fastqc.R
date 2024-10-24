#!/usr/bin/env Rscript
# functions moved

source("../functions/fastqc_parser.R")
source("../functions/data_writer.R")

#' Main FastQC parsing function
main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    
    if (length(args) == 0) {
        log_error("No directory specified")
        stop("Usage: Rscript parse_fastqc.R <directory>")
    }
    
    # Find directory
    directory <- find_data_directory(args[1])
    
    if (is.null(directory)) {
        log_error("Directory not found:", args[1])
        stop("Invalid directory")
    }
    
    # Process FastQC files
    parse_fastqc_files(directory)
}

if (!interactive()) {
    main()
}
