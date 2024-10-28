#!/usr/bin/env Rscript
# R/project_init.R

#' Project Initialization System
#' @description Initialize project environment and dependencies
#' @export

#' Load Project Configuration
#' @param config_path Character Path to configuration file
#' @return List Configuration object
load_project_config <- function(
    config_path = "~/lab_utils/R/config/project_config.R"
) {
    if (!file.exists(config_path)) {
        stop("Configuration file not found: ", config_path)
    }
    source(config_path)
    if (!exists("PROJECT_CONFIG")) {
        stop("PROJECT_CONFIG not defined in configuration file")
    }
    invisible(PROJECT_CONFIG)
}

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

#' Load Core Utilities
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
        # Load configuration first
        config <- load_project_config(config_path)
        
        # Check R version
        check_r_version(
            R_system_version = getRversion(),
            R_project_version = package_version(config$SYSTEM$R_MIN_VERSION)
        )
        
        # Ensure log directory exists
        if (!dir.exists("~/logs")) {
            dir.create("~/logs", recursive = TRUE)
        }
        
        # Load core utilities in order
        load_core_utilities(
            base_path = config$PATHS$R_BASE,
            script_list = c(
                "logging_utils.R",
                "environment_utils.R",
                "validation_utils.R"
            ),
            verbose = verbose
        )
        
        # Initialize logging
        log_file <- initialize_logging(script_name = "project_init")
        log_info("Project initialization started", log_file)
        
        # Log system information
        log_system_info(log_file)
        
        log_info("Project initialization completed successfully", log_file)
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

# Execute if run as script
if (!interactive()) initialize_project_main()
