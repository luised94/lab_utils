#' Package Loading Functions
load_required_packages <- function(packages = CONFIG$REQUIRED_PACKAGES) {
    log_info("Loading required packages")
    
    suppressPackageStartupMessages({
        results <- sapply(packages, require, character.only = TRUE)
    })
    
    if (!all(results)) {
        missing_packages <- packages[!results]
        log_error("Failed to load packages:", paste(missing_packages, collapse = ", "))
        stop("Missing required packages")
    }
    
    return(TRUE)
}

#' Input Validation Functions
validate_input <- function(args) {
    log_info("Validating input arguments")
    
    if (length(args) != 1) {
        log_error("Invalid number of arguments")
        stop("Usage: Rscript script.R <directory_path>")
    }
    
    directory_path <- file.path(CONFIG$PATHS$BASE_DIR, args[1])
    
    if (!dir.exists(directory_path)) {
        log_error("Directory not found:", directory_path)
        stop("Invalid directory")
    }
    
    return(directory_path)
}
