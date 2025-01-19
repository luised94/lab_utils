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

debug_info <- 
    list(
    title = "System Diagnostics",
    "system.hostname" = Sys.info()["nodename"],
    "r.version" = R.version.string,
    "r.platform" = R.version$platform,
    "os.type" = .Platform$OS.type,
    "current.directory" = getwd(),
    "session.timezone" = Sys.timezone(),
    "memory.limit" = memory.limit(),
    "cpu.cores" = parallel::detectCores(),
    "user" = Sys.getenv("USER"),
    "git.branch" = tryCatch({
        system("git rev-parse --abbrev-ref HEAD", intern = TRUE)
    }, error = function(e) "git_not_available")
)

print_debug_info(debug_info)
log_system_diagnostics <- function() {
    return(
        list(
            title = "System Diagnostics",
            "system.hostname" = Sys.info()["nodename"],
            "r.version" = R.version.string,
            "r.platform" = R.version$platform,
            "os.type" = .Platform$OS.type,
            "current.directory" = getwd(),
            "session.timezone" = Sys.timezone(),
            "memory.limit" = memory.limit(),
            "cpu.cores" = parallel::detectCores(),
            "user" = Sys.getenv("USER"),
            "git.branch" = tryCatch({
                system("git rev-parse --abbrev-ref HEAD", intern = TRUE)
            }, error = function(e) "git_not_available")
        )
    )
}
log_session_info <- function() {
    list(
        title = "Session Summary",
        "total_objects" = length(ls(envir = globalenv())),
        "loaded_packages" = names(sessionInfo()$loadedOnly),
        "functions" = ls(envir = globalenv())[
            sapply(ls(envir = globalenv()), 
                   function(x) is.function(get(x, envir = globalenv())))
        ],
        "session_start_time" = Sys.time(),
        "r_version" = R.version.string,
        "r_platform" = R.version$platform
    )
}

# Usage
session_debug_info <- log_session_info()
print_debug_info(session_debug_info)
log_session_info <- function(log_file = NULL) {
    # Capture all objects in global environment
    all_objects <- ls(envir = globalenv())
    
    # Categorize objects
    object_summary <- lapply(all_objects, function(obj) {
        obj_value <- get(obj, envir = globalenv())
        list(
            type = class(obj_value)[1],
            length = length(obj_value),
            memory_size = object.size(obj_value)
        )
    })
    names(object_summary) <- all_objects
    
    # Capture loaded packages
    loaded_packages <- sessionInfo()$loadedOnly
    
    # Capture function names
    functions <- all_objects[sapply(all_objects, function(x) is.function(get(x, envir = globalenv())))]
    
    # Comprehensive session info
    session_info <- list(
        title = "Session Summary",
        "total_objects" = length(all_objects),
        "objects" = object_summary,
        "loaded_packages" = names(loaded_packages),
        "functions" = functions,
        "session_start_time" = Sys.time(),
        "r_session_info" = capture.output(sessionInfo())
    )
    
    # Option to log to file
    if (!is.null(log_file)) {
        writeLines(
            paste(
                "SESSION SUMMARY\n",
                "Total Objects:", length(all_objects), "\n",
                "Loaded Packages:", paste(names(loaded_packages), collapse = ", "), "\n",
                "Functions:", paste(functions, collapse = ", ")
            ), 
            log_file
        )
    }
    
    return(session_info)
}

# Usage at end of script
session_summary <- log_session_info()


print_debug_info(session_summary)


log_system_diagnostics <- function() {
    # Consider adding these useful diagnostics:
    list(
        title = "System Diagnostics",
        # ... existing items ...
            "system.hostname" = Sys.info()["nodename"],
            "r.version" = R.version.string,
            "r.platform" = R.version$platform,
            "os.type" = .Platform$OS.type,
            "current.directory" = getwd(),
            "session.timezone" = Sys.timezone(),
            "memory.limit" = memory.limit(),
            "cpu.cores" = parallel::detectCores(),
            "user" = Sys.getenv("USER"),
            "git.branch" = tryCatch({
                system("git rev-parse --abbrev-ref HEAD", intern = TRUE)
            }, error = function(e) "git_not_available"),
        
        # Additional useful diagnostics:
        "locale" = Sys.getlocale(),
        "temp_dir" = tempdir(),
        "ram_free" = tryCatch({
            as.numeric(system("free -g | grep Mem:", intern = TRUE))
        }, error = function(e) "not_available"),
        "disk_free" = tryCatch({
            system("df -h . | tail -1", intern = TRUE)
        }, error = function(e) "not_available")
    )
}

log_session_info <- function() {  # Remove log_file parameter
    # Get environment objects
    all_objects <- ls(envir = globalenv())
    
    # More efficient object summary
    object_summary <- sapply(all_objects, function(obj) {
        obj_value <- get(obj, envir = globalenv())
        c(
            type = class(obj_value)[1],
            size = format(object.size(obj_value), units = "auto")
        )
    }, simplify = FALSE)
    
    list(
        title = "Session Summary",
        "total_objects" = length(all_objects),
        "objects" = object_summary,
        "loaded_packages" = names(sessionInfo()$loadedOnly),
        "attached_packages" = names(sessionInfo()$otherPkgs),
        "functions" = all_objects[sapply(all_objects, 
            function(x) is.function(get(x, envir = globalenv())))],
        "session_time" = Sys.time()
    )
}
full_diagnostics <- list(
    title = "Complete Diagnostics",
    log_session_info()
)

diagnostics <- log_system_diagnostics()
# Use with print_debug_info
print_debug_info(diagnostics)
session_info <- log_session_info()
print_debug_info(session_info)

