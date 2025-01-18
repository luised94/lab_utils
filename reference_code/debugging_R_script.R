
# Ensure required packages are installed
for (pkg in c("stringi", "jsonlite", "microbenchmark")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, quiet = TRUE)
  }
}

# Load required libraries
suppressPackageStartupMessages({
  library(stringi)
  library(jsonlite)
  library(microbenchmark)
})

# Debug configuration
DEBUG_CONFIG <- new.env()
DEBUG_CONFIG$LEVEL <- 5  # Set to TRACE level
DEBUG_CONFIG$OUTPUT <- "console"  # Default to console output
DEBUG_CONFIG$PERFORMANCE_TRACKING <- TRUE

# Debug levels
DEBUG_LEVELS <- list(
  ERROR = 1,
  WARN = 2,
  INFO = 3,
  DEBUG = 4,
  TRACE = 5
)

# Output handlers
output_handlers <- list(
  console = function(msg) cat(msg, "\n"),
  file = function(msg) {
    cat(msg, "\n", file = "debug_output.log", append = TRUE)
    if (is.list(msg)) {
      msg_json <- toJSON(msg, auto_unbox = TRUE)
      cat(msg_json, "\n", file = "debug_output.json", append = TRUE)
    }
  }
)

# Enhanced print_debug_info function
print_debug_info <- function(info_list, level = DEBUG_LEVELS$DEBUG) {
  if (DEBUG_CONFIG$LEVEL < level) return(invisible(NULL))
  
  start_time <- Sys.time()
  
  # Enhance context
  info_list$system_info <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    pid = Sys.getpid(),
    user = Sys.info()["user"],
    R_version = R.version.string
  )
  
  # Capture call stack
  info_list$call_stack <- capture.output(traceback(2))
  
  # Format output
  output <- format_debug_output(info_list)
  
  # Send to appropriate output handler
  output_handlers[[DEBUG_CONFIG$OUTPUT]](output)
  
  # Performance tracking
  if (DEBUG_CONFIG$PERFORMANCE_TRACKING) {
    end_time <- Sys.time()
    duration <- as.numeric(end_time - start_time, units = "secs")
    performance_log(list(
      function_name = "print_debug_info",
      duration = duration,
      timestamp = as.character(end_time)
    ))
  }
  
  invisible(output)
}

# Helper function to format debug output
format_debug_output <- function(info_list) {
  width <- 80  # Fixed width for consistency
  separator <- stri_dup("=", width)
  
  format_item <- function(key, value, indent = 0) {
    if (is.list(value)) {
      c(sprintf("%s%s:", stri_dup("  ", indent), key),
        unlist(lapply(names(value), function(k) format_item(k, value[[k]], indent + 1))))
    } else {
      sprintf("%s%s: %s", stri_dup("  ", indent), key, as.character(value))
    }
  }
  
  body <- unlist(lapply(names(info_list), function(k) format_item(k, info_list[[k]])))
  
  c(separator, body, separator)
}

# Performance logging
performance_log <- function(log_entry) {
  cat(toJSON(log_entry, auto_unbox = TRUE), "\n", file = "debug_performance.json", append = TRUE)
}

# Sample data processing function
process_sample <- function(sample_id, method) {
  print_debug_info(list(
    title = "Sample Processing Start",
    sample_id = sample_id,
    method = method
  ), level = DEBUG_LEVELS$INFO)
  
  # Simulate processing
  Sys.sleep(0.1)
  
  result <- list(
    sample_id = sample_id,
    processed = TRUE,
    timestamp = Sys.time()
  )
  
  print_debug_info(list(
    title = "Sample Processing Complete",
    result = result
  ), level = DEBUG_LEVELS$DEBUG)
  
  return(result)
}

cat("Starting debug demo...\n")

# Basic usage demonstration
print_debug_info(list(
  title = "Basic Debug Info",
  message = "This is a basic debug message",
  value = 42
))

# Demonstrate different debug levels
print_debug_info(list(title = "Error Level Test"), level = DEBUG_LEVELS$ERROR)
print_debug_info(list(title = "Warn Level Test"), level = DEBUG_LEVELS$WARN)
print_debug_info(list(title = "Info Level Test"), level = DEBUG_LEVELS$INFO)
print_debug_info(list(title = "Debug Level Test"), level = DEBUG_LEVELS$DEBUG)
print_debug_info(list(title = "Trace Level Test"), level = DEBUG_LEVELS$TRACE);

# Demonstrate nested structures
print_debug_info(list(
   title = "Nested Structure Demo",
   config = list(
     database = list(
       host="localhost",
       port=5432,
       credentials=list(
         username="user",
         password="REDACTED"
       )
     ),
     api=list(
       endpoint="https://api.example.com",
       version="v2"
     )
   )
))

# Demonstrate error handling
tryCatch({
   stop("Simulated error")
}, error=function(e) {
   print_debug_info(list(
     title="Error Caught",
     error_message=e$message,
     call=deparse(e$call)
   ), level=DEBUG_LEVELS$ERROR);
})

# Demonstrate file output test directly in console for confirmation.
cat("Testing file output...\n")
DEBUG_CONFIG$OUTPUT <- "file"
print_debug_info(list(
   title="File Output Test",
   message="This message should go to the file"
))
DEBUG_CONFIG$OUTPUT <- "console" # Reset to console for further outputs

# Demonstrate JSON output test directly in console for confirmation.
cat("Testing JSON output...\n")
DEBUG_CONFIG$OUTPUT <- "json"
print_debug_info(list(
   title="JSON Output Test",
   data=list(x=1,y=2,z=list(a=3,b=4))
))
DEBUG_CONFIG$OUTPUT <- "console" # Reset to console for further outputs

# Demonstrate integration with sample processing test directly in console for confirmation.
cat("Processing samples...\n")
samples=c("SAMPLE1","SAMPLE2","SAMPLE3")
methods=c("method_A","method_B","method_C")
results=lapply(seq_along(samples),function(i){
   process_sample(samples[i],methods[i])
})

# Final debug info summary test directly in console for confirmation.
print_debug_info(list(
   title="Processing Summary",
   samples_processed=length(results),
   methods_used=unique(methods)
))

cat("Debug demo completed. Check debug_output.log and debug_output.json for additional output.\n")
