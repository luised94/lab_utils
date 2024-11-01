# tests/test_core.R

#' Test Core Functionality
#' @return Logical TRUE if all tests pass
test_core <- function() {
    # Test logging
    log_file <- initialize_logging("test")
    stopifnot(
        "Log file creation failed" = file.exists(log_file),
        "Log directory incorrect" = grepl(CORE_CONFIG$PATHS$LOGS, log_file)
    )
    
    # Test message writing
    success <- write_log_message("INFO", "Test message", log_file)
    stopifnot(
        "Log writing failed" = success,
        "Message not written" = any(grepl("Test message", readLines(log_file)))
    )
    
    # Test locking
    test_file <- tempfile()
    stopifnot(
        "Lock acquisition failed" = acquire_lock(test_file),
        "Lock release failed" = release_lock(test_file)
    )
    
    return(TRUE)
}

# Run tests
if (!interactive()) {
    if (!test_core()) {
        stop("Core validation failed")
    }
    cat("Core validation successful\n")
}
