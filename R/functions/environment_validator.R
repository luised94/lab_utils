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
#' Environment Validation Functions
validate_environment <- function() {
    log_info("Validating environment")
    
    # Check required environment variables
    check_environment_variables()
    
    # Get user information
    user_info <- get_user_info()
    
    # Validate paths
    validate_paths(user_info)
    
    user_info
}

check_environment_variables <- function() {
    missing_vars <- setdiff(
        CONFIG$ENVIRONMENT$REQUIRED_VARS,
        names(Sys.getenv())
    )
    
    if (length(missing_vars) > 0) {
        log_error("Missing environment variables:", 
                 paste(missing_vars, collapse = ", "))
        log_error("Check bashrc configuration")
        stop("Environment not properly configured")
    }
}

get_user_info <- function() {
    username <- Sys.getenv("WINDOWS_USER")
    
    if (username == "") {
        log_error("WINDOWS_USER not defined")
        stop("Invalid environment configuration")
    }
    
    log_info("Using username:", username)
    
    list(
        username = username,
        home = Sys.getenv("HOME")
    )
}
