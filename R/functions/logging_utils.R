#!/usr/bin/env Rscript

#' Logging Utilities
#' @description Centralized logging system for R scripts
#' @export

#' Initialize logging system
#' @param config List Project configuration
#' @param script_name Character Name of the calling script
#' @return Character Path to log file
initialize_logging <- function(
    config = PROJECT_CONFIG,
    script_name = NULL
) {
    tryCatch({
        # Setup central log directory
        log_dir <- file.path(
            normalizePath("~/logs"),
            format(Sys.Date(), "%Y-%m")
        )
        dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Generate log file path
        log_file <- file.path(
            log_dir,
            paste0(
                format(Sys.Date(), "%Y-%m-%d"),
                "_",
                if (is.null(script_name)) basename(sys.frame(1)$ofile) else script_name,
                ".log"
            )
        )
        
        # Initialize log file with headers
        log_system_info(log_file)
        log_git_info(log_file)
        
        return(log_file)
    }, error = function(e) {
        stop(sprintf("Failed to initialize logging: %s", e$message))
    })
}

#' Log System Information
#' @param log_file Character Path to log file
#' @return None
log_system_info <- function(log_file) {
    log_message(
        level = "INFO",
        message = paste(
            "System Info:",
            "R", R.version.string,
            "on", paste(Sys.info()[c("sysname", "release")], collapse = " ")
        ),
        log_file = log_file
    )
}

#' Log Git Information
#' @param log_file Character Path to log file
#' @return None
log_git_info <- function(log_file) {
    tryCatch({
        git_branch <- system("git rev-parse --abbrev-ref HEAD", intern = TRUE)
        git_hash <- system("git rev-parse HEAD", intern = TRUE)
        log_message(
            level = "INFO",
            message = sprintf("Git: branch=%s commit=%s", git_branch, git_hash),
            log_file = log_file
        )
    }, error = function(e) {
        log_message(
            level = "WARNING",
            message = "Git information unavailable",
            log_file = log_file
        )
    })
}

#' Log Message
#' @param level Character Log level
#' @param message Character Message to log
#' @param log_file Character Path to log file
#' @return None
log_message <- function(
    level,
    message,
    log_file = NULL
) {
    valid_levels <- c("TRACE", "DEBUG", "INFO", "WARNING", "ERROR", "FATAL")
    if (!level %in% valid_levels) {
        stop(sprintf("Invalid log level: %s", level))
    }
    
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    log_entry <- sprintf("[%s] [%s] %s", timestamp, level, message)
    
    # Console output
    cat(log_entry, "\n")
    
    # File output
    if (!is.null(log_file)) {
        write(log_entry, file = log_file, append = TRUE)
    }
}

#' Convenience logging functions
#' @param message Character Message to log
#' @param log_file Character Path to log file
#' @return None
log_trace <- function(message, log_file = NULL) {
    log_message("TRACE", message, log_file)
}

log_debug <- function(message, log_file = NULL) {
    log_message("DEBUG", message, log_file)
}

log_info <- function(message, log_file = NULL) {
    log_message("INFO", message, log_file)
}

log_warning <- function(message, log_file = NULL) {
    log_message("WARNING", message, log_file)
}

log_error <- function(message, log_file = NULL) {
    log_message("ERROR", message, log_file)
}

log_fatal <- function(message, log_file = NULL) {
    log_message("FATAL", message, log_file)
    stop(message)
}
