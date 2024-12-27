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
#!/usr/bin/env Rscript

source("../functions/chromosome_converter.R")
source("../functions/granges_converter.R")

#' Test chromosome conversion functionality
test_chromosome_conversion <- function() {
    log_info("Starting chromosome conversion tests")
    
    # Create test data
    test_gr <- GRanges(
        seqnames = c("1", "X", "chrII", "20", "chrY", "M",
                    "chr10", "III", "chrMT"),
        ranges = IRanges(
            start = sample(1:1000, 9),
            width = 100
        )
    )
    
    # Test conversions
    test_conversions(test_gr)
    
    # Test error handling
    test_error_handling()
    
    log_info("Tests completed")
}

test_conversions <- function(gr) {
    log_info("Testing style conversions")
    
    for (style in names(CONFIG$NAMING$STYLES)) {
        log_info("Testing conversion to:", style)
        result <- convert_granges_style(gr, style, verbose = TRUE)
        print(result)
    }
}

test_error_handling <- function() {
    log_info("Testing error handling")
    
    tryCatch({
        convert_granges_style(data.frame(), "UCSC")
    }, error = function(e) {
        log_info("Expected error caught:", conditionMessage(e))
    })
}

if (!interactive()) {
    test_chromosome_conversion()
}
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
