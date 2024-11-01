#!/usr/bin/env Rscript

#' Logging Utilities
#' @description Centralized logging system for R scripts
#' @export



#' Initialize Logging System
#' @param script_name Character Script identifier
#' @param log_dir Character Optional log directory
#' @return Character Path to log file
initialize_logging <- function(
    script_name = basename(sys.frame(1)$ofile),
    log_dir = file.path(CORE_CONFIG$PATHS$LOGS, format(Sys.Date(), "%Y-%m"))
) {
    tryCatch({

        stopifnot(
            "Context must be character" = is.character(script_name),
            "Log directory must be character" = is.character(log_dir)
        )
        # Create log directory
        dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
        
        # Set log file path
        log_file <- file.path(
            log_dir,
            sprintf(
                "%s_%s.log",
                format(Sys.Date(), "%Y-%m-%d"),
                script_name
            )
        )
        
        # Initialize with atomic write
        if (!file.exists(log_file)) {
            if (acquire_lock(log_file)) {
                writeLines(
                    c(
                        sprintf("=== Log Started: %s ===", format(Sys.time())),
                        sprintf("Script: %s", script_name),
                        sprintf("User: %s", Sys.info()[["user"]]),
                        sprintf("R Version: %s", R.version.string)
                    ),
                    log_file
                )
                release_lock(log_file)
            }
        }
        
        return(log_file)
        
    }, error = function(e) {
        stop(sprintf("Logging initialization failed: %s", e$message))
    })
}


#' Write Log Message
#' @param level Character Log level
#' @param message Character Message to log
#' @param log_file Character Path to log file
#' @return Logical TRUE if successful
write_log_message <- function(
    level,
    message,
    log_file = NULL
) {
    
    # Validate inputs
    stopifnot(
        "Level must be in CORE_CONFIG$LOGGING$LEVELS" = 
            level %in% CORE_CONFIG$LOGGING$LEVELS,
        "Message must be character" = is.character(message),
        "Log file must be character" = is.character(log_file) || is.null(log_file)
    )
    # Format message
    entry <- sprintf(
        "[%s] [%s] %s",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        level,
        message
    )
    
    # Console output
    cat(entry, "\n", file = stderr())
    
    # File output with lock
    if (!is.null(log_file)) {
        if (acquire_lock(log_file)) {
            write(entry, file = log_file, append = TRUE)
            release_lock(log_file)
            return(TRUE)
        }
        return(FALSE)
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
write_log_trace <- function(message, log_file = NULL) {
    write_log_message("TRACE", message, log_file)
}

write_log_debug <- function(message, log_file = NULL) {
    write_log_message("DEBUG", message, log_file)
}

write_log_info <- function(message, log_file = NULL) {
    write_log_message("INFO", message, log_file)
}

write_log_warning <- function(message, log_file = NULL) {
    write_log_message("WARNING", message, log_file)
}

write_log_error <- function(message, log_file = NULL) {
    write_log_message("ERROR", message, log_file)
}

write_log_fatal <- function(message, log_file = NULL) {
    write_log_message("FATAL", message, log_file)
    stop(message)
}
