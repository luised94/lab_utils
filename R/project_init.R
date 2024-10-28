#!/usr/bin/env Rscript

#' Project Initialization Main Function
#' @param config_path Character. Path to configuration file
#' @param force Logical. Force reinitialization if project exists
#' @param verbose Logical. Enable detailed logging
#' @return Logical. TRUE if initialization successful
#' @export
initialize_project_main <- function(
    config_path = "R/config/project_config.R",
    force = FALSE,
    verbose = TRUE
) {
    # Setup error handling with context
    options(error = function() {
        message("Error in project initialization")
        if (interactive()) recover()
    })

    tryCatch({

        # Validate R version
        check_r_version(
            R_system_version = getRversion(),
            R_project_version = package_version(PROJECT_CONFIG$SYSTEM$R_MIN_VERSION)
        )
        
        # Load and validate configuration
        config <- load_project_config(config_path)
        
        # Initialize directory structure
        setup_project_structure(
            root_dir = config$PATHS$ROOT,
            force = force
        )
        
        # Initialize renv
        setup_renv_environment(
            project_root = config$PATHS$ROOT,
            settings = config$ENVIRONMENT$RENV_SETTINGS
        )
        
        # Load core utilities
        load_core_utilities(
            base_path = config$PATHS$R_BASE,
            scripts = config$LOAD_SEQUENCE$CRITICAL,
            verbose = verbose
        )
        
        message("Project initialization completed successfully")
        invisible(TRUE)
        
    }, error = function(e) {
        message("Critical initialization error: ", e$message)
        if (interactive()) {
            recover()
        } else {
            quit(status = 1)
        }
    })
}

#' Check R Version Compatibility
#' @return Logical. TRUE if version check passes
check_r_version <- function(R_system_version, R_project_version) {
    current_version <- R_system_version
    required_version <- R_project_version
    
    if (current_version < required_version) {
        stop(sprintf(
            "R version %s required, current version is %s",
            required_version,
            current_version
        ))
    }
    invisible(TRUE)
}

#' Setup Project Directory Structure
#' @param root_dir Character. Project root directory
#' @param force Logical. Force recreation of directories
#' @return Logical. TRUE if setup successful
setup_project_structure <- function(root_dir, force = FALSE) {
    required_dirs <- c("R", "logs", "data", "output")
    
    for (dir in required_dirs) {
        dir_path <- file.path(root_dir, dir)
        if (!dir.exists(dir_path)) {
            dir.create(dir_path, recursive = TRUE)
        }
    }
    invisible(TRUE)
}

#' Initialize renv Environment
#' @param project_root Character. Project root directory
#' @param settings List. renv configuration settings
#' @return Logical. TRUE if initialization successful
setup_renv_environment <- function(project_root, settings) {
    if (!requireNamespace("renv", quietly = TRUE)) {
        install.packages("renv")
    }
    
    renv::init(
        project = project_root,
        force = settings$AUTO_ACTIVATE
    )
    invisible(TRUE)
}

# Execute if run as script
if (!interactive()) {
    initialize_project_main()
}
