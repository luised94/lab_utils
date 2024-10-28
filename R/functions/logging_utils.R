#' Advanced Logging Functions for R Scripts
#'
#' @title Advanced logging utility functions
#' @description A set of functions for consistent logging across R scripts
#' @author Your Name
#' @date 2024-10-18

library(tools)
source("~/lab_utils/R/init.R")

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

scp_command_reminder <- function(verbose = TRUE) {

    cat("Ensure you are connected to the luria mit network.\n")
    cat("Run the following command to rsync the created directory.\n")
    server_path <- "luised94@luria.mit.edu:~/data/"
    cat("scp -r from_dir user@server:to_dir\n")
    cat(sprintf("scp -r \"%s\" \"%s\"", experiment_dir, server_path), "\n")
    cat("After running the scp command, login to cluster and \n download the data from BMC (see 001_downloadDataFromBMC.sh )\n")
    print("Script complete.")
}







