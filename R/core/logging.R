#!/usr/bin/env Rscript

#' Logging Utilities
#' @description Centralized logging system for R scripts
#' @export


#' Initialize Logging System
#' @param context Character Logging context identifier
#' @param log_dir Character Optional custom log directory
#' @return Character Path to log file
#' @export
initialize_logging <- function(
    context = "main",
    log_dir = CORE_CONFIG$LOGGING$DEFAULT_ROOT
) {
    # Validate inputs
    stopifnot(
        "Context must be character" = is.character(context),
        "Log directory must be character" = is.character(log_dir)
    )
    
    # Create log directory
    log_path <- file.path(
        log_dir,
        format(Sys.time(), "%Y-%m")
    )
    dir.create(log_path, recursive = TRUE, showWarnings = FALSE)
    
    # Set log file
    log_file <- file.path(
        log_path,
        sprintf("%s_%s.log", 
                format(Sys.time(), "%Y-%m-%d"),
                context)
    )
    
    # Initialize if new
    if (!file.exists(log_file)) {
        cat(sprintf(
            "=== Log Initialized: %s ===\n",
            format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        ), file = log_file)
    }
    
    return(log_file)
}

#' Write Log Message
#' @param level Character Log level
#' @param message Character Message to log
#' @param log_file Character Path to log file
#' @param context Character Optional context
#' @return Logical TRUE if successful
#' @export
write_log_message <- function(
    level = "INFO",
    message,
    log_file,
    context = "main"
) {
    # Validate inputs
    stopifnot(
        "Level must be in CORE_CONFIG$LOGGING$LEVELS" = 
            level %in% CORE_CONFIG$LOGGING$LEVELS,
        "Message must be character" = is.character(message),
        "Log file must be character" = is.character(log_file)
    )
    
    # Format message
    log_entry <- sprintf(
        CORE_CONFIG$LOGGING$FORMAT,
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        level,
        context,
        message
    )
    
    # Write with lock
    acquire_lock(log_file) && {
        cat(log_entry, "\n", file = log_file, append = TRUE)
        release_lock(log_file)
    }
    
    return(TRUE)
}

#' Check Log Status
#' @param log_file Character Path to log file
#' @return List Log status information
check_log_status <- function(log_file) {
    list(
        exists = file.exists(log_file),
        size = if (file.exists(log_file)) file.info(log_file)$size else 0,
        last_modified = if (file.exists(log_file)) 
            format(file.info(log_file)$mtime, "%Y-%m-%d %H:%M:%S") else NA
    )
}


#' Log System Information
#' @param log_file Character Path to log file
#' @return None
log_system_info <- function(log_file) {
    write_log_message(
        level = "INFO",
        message = paste(
            "System Info:",
            "R", R.version.string,
            "on", paste(Sys.info()[c("sysname", "release")], collapse = " "),
            "Platform:", sessionInfo()$platform
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
        write_log_message(
            level = "INFO",
            message = sprintf("Git: branch=%s commit=%s", git_branch, git_hash),
            log_file = log_file
        )
    }, error = function(e) {
        write_log_message(
            level = "WARNING",
            message = "Git information unavailable",
            log_file = log_file
        )
    })
}


#' Convenience logging functions
#' @param message Character Message to log
#' @param log_file Character Path to log file
#' @return None
log_trace <- function(message, log_file = NULL) {
    write_log_message("TRACE", message, log_file)
}

log_debug <- function(message, log_file = NULL) {
    write_log_message("DEBUG", message, log_file)
}

log_info <- function(message, log_file = NULL) {
    write_log_message("INFO", message, log_file)
}

log_warning <- function(message, log_file = NULL) {
    write_log_message("WARNING", message, log_file)
}

log_error <- function(message, log_file = NULL) {
    write_log_message("ERROR", message, log_file)
}

log_fatal <- function(message, log_file = NULL) {
    write_log_message("FATAL", message, log_file)
    stop(message)
}
