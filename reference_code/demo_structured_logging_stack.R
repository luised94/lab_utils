#!/usr/bin/env Rscript

# Minimal demonstration of a structured logging approach in R

# 1. Library Imports
library(futile.logger)  # Provides the flog.* logging functions
library(stringi)        # For string manipulation (e.g., stri_dup)
library(crayon)         # For colored output, if desired

# 2. Configure Logging
flog.threshold(DEBUG)  # Ensure all log levels (DEBUG, INFO, etc.) are shown

# 3. Structured Logging Function
structured_log_info <- function(message, step = NULL) {
  # Retrieve the call stack
  calls <- sys.calls()
  # The second-to-last call typically corresponds to the function that invoked
  # structured_log_info (i.e., user-level code).
  caller_call <- calls[[length(calls) - 1]]
  function_context <- as.character(caller_call)
  
  # Build a timestamp
  timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")
  
  # Build a logfmt-like message for easy parsing with grep, sed, awk, etc.
  structured_message <- paste0(
    "level=INFO ",
    "time=", timestamp, " ",
    "function=", function_context, " ",
    if (!is.null(step)) paste0("step=", step, " ") else "",
    'message="', message, '"'
  )
  
  # Log the message using futile.logger
  flog.info(structured_message)
}

# 4. Example Utility Functions
get_width <- function() {
  width <- Sys.getenv("COLUMNS")
  if (width != "") {
    as.integer(width)
  } else {
    tryCatch({
      as.integer(system("tput cols", intern = TRUE))
    }, error = function(e) {
      80  # Default width if tput is not available
    })
  }
}

create_separator <- function(char = "=", width = get_width()) {
  stri_dup(char, width)
}

# 5. Example Debug Info List (Not used by a printing function here; just a placeholder)
debug_info <- list(
  "title" = "Debug Information",
  "Scaling Mode" = "Some Mode",
  "Plot Generation Details" = NULL,
  ".Experiment" = "Exp001",
  ".Chromosome" = "Chr1",
  ".Sample Count" = 10,
  ".Timestamp" = Sys.time(),
  ".Normalization" = "Method1",
  "Output Configuration" = NULL,
  ".Plot Directory" = "/path/to/output",
  ".Filename" = "plot.png",
  ".Full Path" = "/path/to/output/plot.png"
)

# 6. Example Functions Demonstrating Stack Context
#    Each function logs a line on entry or exit using structured_log_info.

f1 <- function() {
  structured_log_info("Entered f1", step = "init")
  f2()
  structured_log_info("Exiting f1", step = "done")
}

f2 <- function() {
  structured_log_info("Entered f2", step = "processing")
  # Pretend we're doing some work here
  structured_log_info("Finished processing in f2", step = "completion")
}

# 7. Demonstration of Use
structured_log_info("Starting demonstration script")
structured_log_info("Now calling f1", step = "main")

f1()

structured_log_info("Demonstration script complete", step = "main")
