# tests/test_core.R

#' Test Core Functionality
#' @return Logical TRUE if all tests pass
test_core <- function() {
    # Test logging
    log_file <- initialize_logging(
        context = "test",
        log_dir = tempdir()
    )
    
    stopifnot(
        "Log file creation failed" = file.exists(log_file),
        "Log message failed" = write_log_message(
            level = "INFO",
            message = "Test message",
            log_file = log_file
        )
    )
    
    # Test locking
    test_file <- tempfile()
    stopifnot(
        "Lock acquisition failed" = acquire_lock(test_file),
        "Lock release failed" = release_lock(test_file)
    )
    
    return(TRUE)
}

# Quick test
if (!interactive()) {
    if (!test_core()) {
        stop("Core validation failed")
    }
    cat("Core validation successful\n")
}
