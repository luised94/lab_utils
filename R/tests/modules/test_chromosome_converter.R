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
