#' Environment Setup Functions
initialize_project_environment <- function(config = CONFIG$ENVIRONMENT) {
    log_info("Initializing project environment")
    
    # Validate environment
    validate_environment(config)
    
    # Setup paths
    setup_project_paths(config)
    
    # Initialize renv if needed
    initialize_renv_if_needed()
    
    # Load required packages
    load_required_packages()
}

setup_project_paths <- function(config) {
    log_info("Setting up project paths")
    
    # Ensure R library path
    dir.create(
        config$PATHS$USER_LIB,
        recursive = TRUE,
        showWarnings = FALSE
    )
    
    # Set library paths
    .libPaths(c(
        config$PATHS$USER_LIB,
        .libPaths()
    ))
}

initialize_renv_if_needed <- function() {
    log_info("Checking renv status")
    
    if (!requireNamespace("renv", quietly = TRUE)) {
        log_info("Initializing renv")
        install.packages("renv")
        renv::init()
    }
}
