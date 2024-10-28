#!/usr/bin/env Rscript
# R/functions/environment_utils.R

#' Environment Management Utilities
#' @description Centralized environment management system
#' @export

#' Validate System Environment
#' @param config List Configuration settings
#' @param min_memory Numeric Minimum required memory in GB
#' @return List Environment status
validate_system_environment <- function(
    config = PROJECT_CONFIG,
    min_memory = 4
) {
    tryCatch({
        # Check R version
        check_r_version(
            R_system_version = getRversion(),
            R_project_version = package_version(config$SYSTEM$R_MIN_VERSION)
        )
        
        # Check memory
        check_system_memory(
            min_memory = min_memory,
            current_memory = memory.limit()
        )
        
        # Check environment variables
        check_environment_variables(
            required_vars = config$ENVIRONMENT$REQUIRED_VARS,
            current_vars = names(Sys.getenv())
        )
        
        # Log system info
        log_system_info()
        
        invisible(TRUE)
    }, error = function(e) {
        log_error("Environment validation failed:", e$message)
        stop(e)
    })
}

#' Check System Memory
#' @param min_memory Numeric Required memory in GB
#' @param current_memory Numeric Available memory in MB
#' @return Logical TRUE if check passes
check_system_memory <- function(
    min_memory,
    current_memory
) {
    if (current_memory < min_memory * 1024) {
        log_warning(sprintf(
            "Low memory: %d MB available, %d GB required",
            current_memory,
            min_memory
        ))
    }
    invisible(TRUE)
}

#' Check Environment Variables
#' @param required_vars Character vector Required variables
#' @param current_vars Character vector Current variables
#' @return Logical TRUE if check passes
check_environment_variables <- function(
    required_vars,
    current_vars
) {
    missing_vars <- setdiff(required_vars, current_vars)
    if (length(missing_vars) > 0) {
        stop(sprintf(
            "Missing environment variables: %s",
            paste(missing_vars, collapse = ", ")
        ))
    }
    invisible(TRUE)
}

#' Log System Information
#' @return None
log_system_info <- function() {
    info <- sessionInfo()
    log_info(sprintf("R Version: %s", info$R.version$version.string))
    log_info(sprintf("Platform: %s", info$platform))
    log_info(sprintf("Working Directory: %s", getwd()))
}

#' Setup Project Paths
#' @param config List Configuration settings
#' @return List Path information
setup_project_paths <- function(
    config = PROJECT_CONFIG
) {
    tryCatch({
        # Validate base paths
        for (path in c(config$PATHS$ROOT, "~/logs")) {
            if (!dir.exists(path)) {
                dir.create(path, recursive = TRUE)
            }
        }
        
        # Return path information
        list(
            root = normalizePath(config$PATHS$ROOT),
            logs = normalizePath("~/logs"),
            data = if (dir.exists("~/data")) normalizePath("~/data") else NULL
        )
    }, error = function(e) {
        log_error("Path setup failed:", e$message)
        stop(e)
    })
}

#' Clean R Environment
#' @param keep Character vector Names to keep
#' @return None
clean_environment <- function(
    keep = c("PROJECT_CONFIG")
) {
    # Detach packages
    attached <- paste0("package:", names(sessionInfo()$otherPkgs))
    for (pkg in attached) {
        try(detach(pkg, character.only = TRUE, unload = TRUE), silent = TRUE)
    }
    
    # Clear workspace except kept objects
    rm(list = setdiff(ls(all.names = TRUE), keep), envir = .GlobalEnv)
    
    # Force garbage collection
    gc()
    
    invisible(TRUE)
}
