#' Advanced Logging Functions for R Scripts
#'
#' @title Advanced logging utility functions
#' @description A set of functions for consistent logging across R scripts
#' @author Your Name
#' @date 2024-10-18

library(tools)

#' @title Get script name
#' @description Extract the full path of the current script
#' @return Character string of the script path
get_script_name <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile))
  } else {
    return("interactive")
  }
}

#' @title Get script directory
#' @description Extract the directory of the current script
#' @return Character string of the script directory
get_script_dir <- function() {
  script_path <- get_script_name()
  if (script_path != "interactive") {
    return(dirname(script_path))
  } else {
    return(getwd())
  }
}

#' @title Get script basename
#' @description Extract the basename of the current script without extension
#' @return Character string of the script basename
get_script_basename <- function() {
  script_path <- get_script_name()
  if (script_path != "interactive") {
    return(tools::file_path_sans_ext(basename(script_path)))
  } else {
    return("interactive")
  }
}

#' @title Initialize logging
#' @description Set up logging for a script
#' @param log_file Optional. Path to the log file. If NULL, a default is created.
#' @return The path to the log file
init_logging <- function(log_file = NULL) {
  if (is.null(log_file)) {
    script_name <- get_script_basename()
    log_dir <- file.path(get_script_dir(), "logs", format(Sys.Date(), "%Y-%m"))
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    log_file <- file.path(log_dir, paste0(format(Sys.Date(), "%Y-%m-%d"), "_", script_name, ".log"))
  }
  
  log_system_info(log_file)
  log_git_info(log_file)
  
  return(log_file)
}

#' @title Log system information
#' @description Log R version and system information
#' @param log_file Path to the log file
#' @return None
log_system_info <- function(log_file) {
  log_message("INFO", paste("R version:", R.version.string), log_file)
  log_message("INFO", paste("System:", Sys.info()["sysname"], Sys.info()["release"]), log_file)
}

#' @title Log git information
#' @description Log git branch and commit hash
#' @param log_file Path to the log file
#' @return None
log_git_info <- function(log_file) {
  tryCatch({
    git_branch <- system("git rev-parse --abbrev-ref HEAD", intern = TRUE)
    git_hash <- system("git rev-parse HEAD", intern = TRUE)
    log_message("INFO", paste("Git branch:", git_branch), log_file)
    log_message("INFO", paste("Git commit:", git_hash), log_file)
  }, error = function(e) {
    log_message("WARNING", "Failed to retrieve git information", log_file)
  })
}

#' @title Log a message
#' @description Log a message with timestamp and level
#' @param level The log level (e.g., "INFO", "WARNING", "ERROR")
#' @param message The message to log
#' @param log_file Path to the log file
#' @return None
log_message <- function(level, message, log_file = NULL) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_entry <- sprintf("[%s] [%s] %s", timestamp, level, message)
  
  cat(log_entry, "\n")
  
  if (!is.null(log_file)) {
    write(log_entry, file = log_file, append = TRUE)
  }
}

#' @title Log an info message
#' @description Log a message at the INFO level
#' @param message The message to log
#' @param log_file Path to the log file
#' @return None
log_info <- function(message, log_file = NULL) {
  log_message("INFO", message, log_file)
}

#' @title Log a warning message
#' @description Log a message at the WARNING level
#' @param message The message to log
#' @param log_file Path to the log file
#' @return None
log_warning <- function(message, log_file = NULL) {
  log_message("WARNING", message, log_file)
}

#' @title Log an error message
#' @description Log a message at the ERROR level
#' @param message The message to log
#' @param log_file Path to the log file
#' @return None
log_error <- function(message, log_file = NULL) {
  log_message("ERROR", message, log_file)
}
