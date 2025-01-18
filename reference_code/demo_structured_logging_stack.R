#!/usr/bin/env Rscript

# Minimal demonstration of a structured logging approach in R

# 1. Library Imports
library(futile.logger)  # Provides the flog.* logging functions
library(stringi)        # For string manipulation (e.g., stri_dup)
library(crayon)         # For colored output, if desired

# 2. Configure Logging
flog.threshold(DEBUG)  # Ensure all log levels (DEBUG, INFO, etc.) are shown

# 3. Structured Logging Function

library(futile.logger)
structured_log_info <- function(
  message,
  step = NULL,
  skip_functions = c("eval", "envir", "source", "withVisible", "structured_log_info"),
  log_level = "INFO"
) {
  # --------------------------------------------------------------------------
  # 1. Input Validation
  # --------------------------------------------------------------------------
  if (!is.character(message) || length(message) != 1) {
    stop("`message` must be a single character string.")
  }
  if (!is.null(step) && (!is.character(step) || length(step) != 1)) {
    stop("`step`, if provided, must be a single character string.")
  }
  if (!is.character(skip_functions)) {
    stop("`skip_functions` must be a character vector of function names to ignore.")
  }
  if (!is.character(log_level) || length(log_level) != 1) {
    stop("`log_level` must be a single character string (e.g., \"INFO\", \"DEBUG\").")
  }
  
  # --------------------------------------------------------------------------
  # 2. Capture & Prune the Call Stack
  # --------------------------------------------------------------------------
  calls <- sys.calls()
  
  # Pull out the function names from each call
  call_names <- sapply(calls, function(a_call) {
    func_part <- a_call[[1]]
    as.character(func_part)
  })
  
  # Reverse the vector so the most recent call is first
  reversed_names <- rev(call_names)
  
  # Identify the first user-defined call by skipping known base / wrapper calls
  idx_user_call <- which(!reversed_names %in% skip_functions)[1]
  function_context <- "unknown_function"
  if (!is.na(idx_user_call)) {
    # Convert index in reversed array to index in original array
    actual_index <- length(reversed_names) - idx_user_call + 1
    function_context <- call_names[actual_index]
  }
  
  # --------------------------------------------------------------------------
  # 3. Construct the Log Message
  # --------------------------------------------------------------------------
  timestamp <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")
  
  # Build a logfmt-like message for easy grep/sed/awk usage
  # If you prefer JSON, you could build a JSON string here.
  structured_message <- paste0(
    "level=", log_level, " ",
    "time=", timestamp, " ",
    "function=", function_context, " ",
    if (!is.null(step)) paste0("step=", step, " ") else "",
    "message=\"", message, "\""
  )
  
  # --------------------------------------------------------------------------
  # 4. Log Using the Appropriate Level
  # --------------------------------------------------------------------------
  # We can dispatch based on `log_level` to use flog.debug(), flog.info(), etc.,
  # or we can always use flog.info and rely on the textual level.  
  # Here we match the actual futile.logger function to the parameter.
  
  switch(
    toupper(log_level),
    "DEBUG" = flog.debug(structured_message),
    "INFO"  = flog.info(structured_message),
    "WARN"  = flog.warn(structured_message),
    "ERROR" = flog.error(structured_message),
    # Default if unrecognized log level
    flog.info(structured_message)
  )
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
