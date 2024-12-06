#!/usr/bin/env Rscript
# functions moved

source("../functions/data_downloader.R")
source("../functions/feature_processor.R")

#' Main feature download function
main <- function() {
    log_info("Starting feature download")
    
    # Setup
    timestamp <- format(Sys.time(), "%y%m%d")
    base_dir <- file.path(Sys.getenv("HOME"), "data", "feature_files")
    
    ensure_directory(base_dir)
    
    # Download and process data
    results <- download_feature_data(
        CONFIG$SOURCES,
        base_dir,
        timestamp
    )
    
    # Validate results
    validate_results(results)
    
    log_info("Feature download complete")
    log_info("Cleanup hint: rm *_eaton_acs.bed if not needed")
}

if (!interactive()) {
    main()
}
