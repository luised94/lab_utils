#!/usr/bin/env Rscript
# functions moved

source("../functions/feature_processor.R")
source("../functions/data_converter.R")

#' Main feature processing function
main <- function() {
    log_info("Starting feature processing")
    
    # Process features
    results <- process_feature_files()
    
    if (is.null(results)) {
        log_error("Feature processing failed")
        quit(status = 1)
    }
    
    # Output transfer command
    log_info("Processing complete")
    log_info("To transfer files:")
    log_info("scp -r user@server:from_dir to_dir")
}

if (!interactive()) {
    main()
}
