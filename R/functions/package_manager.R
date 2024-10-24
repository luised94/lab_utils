#' Package Management Functions
load_required_packages <- function(packages = CONFIG$PACKAGES$REQUIRED) {
    log_info("Loading required packages")
    
    suppressPackageStartupMessages({
        results <- sapply(packages, require, character.only = TRUE)
    })
    
    missing <- packages[!results]
    if (length(missing) > 0) {
        log_error("Missing packages:", paste(missing, collapse = ", "))
        stop("Required packages not available")
    }
    
    return(TRUE)
}
