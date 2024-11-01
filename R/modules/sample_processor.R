#' Sample Table Processing Functions
process_control_factors <- function(sample_table) {
    log_info("Processing control factors")
    
    cf_cols <- grep(CONFIG$PATTERNS$CONTROL_FACTOR_PREFIX, names(sample_table), 
                   value = TRUE)
    
    if (length(cf_cols) == 0) {
        log_error("No control factor columns found")
        stop("Invalid sample table format")
    }
    
    control_factors <- lapply(sample_table[cf_cols], function(x) {
        strsplit(x[1], ",")[[1]]
    })
    
    names(control_factors) <- sub(CONFIG$PATTERNS$CONTROL_FACTOR_PREFIX, "", 
                                cf_cols)
    
    sample_table[cf_cols] <- NULL
    attr(sample_table, "control_factors") <- control_factors
    
    return(sample_table)
}

#' Control Sample Matching Functions
determine_matching_control <- function(sample_row, sample_table, factors_to_match) {
    log_info("Finding matching control sample")
    
    comparison_row <- sample_row[factors_to_match]
    matches <- apply(sample_table[, factors_to_match], 1, function(row) {
        all(row == comparison_row)
    })
    
    control_indices <- which(matches & sample_table$antibody == "Input")
    
    if (length(control_indices) == 0) {
        log_warning("No matching control found")
        return(1)
    }
    
    if (length(control_indices) > 1) {
        log_warning("Multiple controls found, using first")
    }
    
    return(control_indices[1])
}
