source("~/lab_utils/core_scripts/functions_for_logging.R")
source("~/lab_utils/core_scripts/template_bmc_config.R")

# Set up the logger to log to a file
#flog.appender(appender.file("debug_log.txt"))  # Change filename as needed
flog.threshold(DEBUG)  # Set minimum log level to DEBUG
# Example usage
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
  "..Full Path" = "/path/to/output/plot.png"
)
print_debug_info(debug_info)

# Set up the logger
#flog.appender(appender.file("logfile.txt"))  # Log to a file
flog.threshold(DEBUG)  # Set minimum log level to INFO

# Example of logging at different levels
flog.debug("This is a debug message.")
flog.info("This is an info message.")
flog.warn("This is a warning message.")
flog.error("This is an error message.")
print_config_settings(RUNTIME_CONFIG)
print_config_settings(GENOME_TRACK_CONFIG)
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

print_debug_info(debug_info)

diagnostics <- log_system_diagnostics()
# Use with print_debug_info
print_debug_info(diagnostics)
session_info <- log_session_info()
print_debug_info(session_info)
combined_diagnostics <- modifyList(
    log_system_diagnostics(),
    log_session_info()
)
print_debug_info(combined_diagnostics)
