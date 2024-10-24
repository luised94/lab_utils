#' Environment Management Functions
check_environment <- function() {
    log_info("Checking environment")
    
    # Check required variables
    missing_vars <- setdiff(
        CONFIG$ENVIRONMENT$REQUIRED_VARS,
        names(Sys.getenv())
    )
    
    if (length(missing_vars) > 0) {
        log_warning("Missing environment variables:",
                   paste(missing_vars, collapse = ", "))
    }
    
    # Set basic options
    setup_r_options()
}

setup_r_options <- function() {
    if (interactive()) {
        options(
            width = 120,
            prompt = "R> ",
            continue = "... ",
            digits = 4,
            show.signif.stars = FALSE
        )
    }
}

get_script_dir <- function() {
    # Get the directory of the currently executing script
    args <- commandArgs(trailingOnly = FALSE)
    script_path <- args[grep("--file=", args)]
    
    if (length(script_path) > 0) {
        # Remove "--file=" prefix
        script_path <- substring(script_path, 8)
        return(dirname(script_path))
    }
    
    # Fallback to current directory
    getwd()
}
