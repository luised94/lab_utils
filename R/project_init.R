#!/usr/bin/env Rscript

#' Check R Version Compatibility
#' @param R_system_version package_version Current R version
#' @param R_project_version package_version Required R version
#' @return Logical TRUE if version check passes
check_r_version <- function(
    R_system_version,
    R_project_version
) {
    if (R_system_version < R_project_version) {
        stop(sprintf(
            "R version %s required, current version is %s",
            R_project_version,
            R_system_version
        ))
    }
    invisible(TRUE)
}

#' Load Project Configuration
#' @param config_path Character Path to configuration file
#' @return List Configuration object
load_project_config <- function(config_path) {
    if (!file.exists(config_path)) {
        stop("Configuration file not found: ", config_path)
    }
    source(config_path)
    return(PROJECT_CONFIG)
}

#' Load Core Utility Functions
#' @param base_path Character Base path for utilities
#' @param script_list Character vector List of scripts to load
#' @param verbose Logical Enable detailed logging
#' @return Logical TRUE if loading successful
load_core_utilities <- function(
    base_path,
    script_list,
    verbose = TRUE
) {
    for (script in script_list) {
        script_path <- file.path(base_path, "functions", script)
        if (!file.exists(script_path)) {
            stop("Required utility not found: ", script_path)
        }
        if (verbose) message("Loading: ", script_path)
        source(script_path)
    }
    invisible(TRUE)
}

#' Project Initialization Main Function
#' @param config_path Character Path to configuration file
#' @param verbose Logical Enable detailed logging
#' @return Logical TRUE if initialization successful
#' @export
initialize_project_main <- function(
    config_path = "~/lab_utils/R/config/project_config.R",
    verbose = TRUE
) {
    tryCatch({
        # Create logs directory if it doesn't exist
        if (!dir.exists("~/logs")) {
            dir.create("~/logs", recursive = TRUE)
        }
        
        # Load configuration
        config <- load_project_config(config_path)
        
        # Validate R version
        check_r_version(
            R_system_version = getRversion(),
            R_project_version = package_version(config$SYSTEM$R_MIN_VERSION)
        )
        
        # Load core utilities
        load_core_utilities(
            base_path = config$PATHS$R_BASE,
            script_list = config$LOAD_SEQUENCE$CRITICAL,
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

if (!interactive()) initialize_project_main()


#
#' Initialize renv Environment
#' @param project_root Character. Project root directory
#' @param settings List. renv configuration settings
#' @return Logical. TRUE if initialization successful
setup_renv_environment <- function(project_root, settings) {
  if (!requireNamespace("renv", quietly = TRUE)) {
    message("Installing renv package...")
    install.packages("renv")
  }
  
  renv_lockfile <- file.path(project_root, "renv.lock")
  
  if (file.exists(renv_lockfile)) {
    message("Activating existing renv environment...")
    renv::activate(project = project_root)
  } else {
    message("Initializing new renv environment...")
    renv::init(project = project_root, force = settings$AUTO_ACTIVATE)
  }
  
  message("renv setup complete.")
  invisible(TRUE)
}











