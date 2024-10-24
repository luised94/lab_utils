#' Environment Validation Functions
validate_environment <- function() {
    log_info("Validating R environment")
    
    # Check R version
    validate_r_version()
    
    # Check available memory
    validate_memory()
    
    # Check package installation
    validate_packages()
    
    # Print session info
    log_session_info()
}

validate_r_version <- function(min_version = "4.0.0") {
    current <- getRversion()
    if (current < min_version) {
        log_warning("R version might be too old:",
                   current, "<", min_version)
    }
}

validate_memory <- function(min_gb = 4) {
    mem <- memory.limit()
    if (mem < min_gb * 1024) {
        log_warning("Available memory might be too low:",
                   mem, "MB")
    }
}

validate_packages <- function() {
    installed <- installed.packages()
    missing <- setdiff(
        unlist(CONFIG$PACKAGES),
        rownames(installed)
    )
    
    if (length(missing) > 0) {
        log_warning("Missing packages:",
                   paste(missing, collapse = ", "))
    }
}

log_session_info <- function() {
    log_info("Session Information:")
    print(sessionInfo())
}
