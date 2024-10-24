#!/usr/bin/env Rscript

source("../functions/input_validator.R")

#' Main plotting function
main <- function() {
    # Get command line arguments
    args <- commandArgs(trailingOnly = TRUE)
    
    # Validate inputs
    tryCatch({
        validate_script_args(
            directory = args[1],
            chromosome = if(length(args) > 1) as.integer(args[2]) else NULL
        )
    }, error = function(e) {
        log_error("Validation failed:", e$message)
        quit(status = 1)
    })
    
    # Continue with script execution...
}

if (!interactive()) {
    main()
}
