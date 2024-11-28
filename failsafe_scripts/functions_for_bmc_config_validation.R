
validate_category_values <- function(
    categories,     # list: category name -> character vector mapping
    stop_on_error = TRUE,  # logical: whether to stop or return validation status
    verbose = FALSE    # logical: whether to print validation steps
) {
    tryCatch({
        if (verbose) cat("[VALIDATING] Checking category values...\n")
        
        sapply(names(categories), function(category_name) {
            values <- categories[[category_name]]
            
            if (verbose) {
                cat(sprintf("  Checking category '%s':\n", category_name))
                cat(sprintf("    Values: %s\n", paste(values, collapse = ", ")))
            }
            
            # Check if values are character
            if (!is.character(values)) {
                if (verbose) cat(sprintf("    [FAIL] Not character vector\n"))
                stop(sprintf(
                    "Category '%s' values must be character vectors",
                    category_name
                ))
            }
            if (verbose) cat("    [PASS] Character vector check\n")
            
            # Check for duplicates
            if (any(duplicated(values))) {
                if (verbose) cat("    [FAIL] Contains duplicates\n")
                stop(sprintf(
                    "Category '%s' values must be unique",
                    category_name
                ))
            }
            if (verbose) cat("    [PASS] Uniqueness check\n")
        })
        
        if (verbose) cat("[PASS] All category values valid\n")
        return(TRUE)
    }, error = function(e) {
        if (verbose) cat(sprintf("[FAIL] Category validation: %s\n", e$message))
        if (stop_on_error) stop(e) else return(FALSE)
    })
}

validate_column_references <- function(
    categories,        # list: valid column names
    comparisons,       # list: named expressions
    control_factors,   # list: factor definitions
    conditions,        # list: experimental conditions
    stop_on_error = TRUE,
    verbose = FALSE
) {
    valid_columns <- names(categories)
    
    if (verbose) {
        cat("[VALIDATING] Checking column references...\n")
        cat("  Valid columns:", paste(valid_columns, collapse = ", "), "\n")
    }
    
    # Helper function for expression validation
    check_expr_vars <- function(expr, context_name) {
        if (verbose) cat(sprintf("  Checking %s\n", context_name))
        
        invalid_cols <- setdiff(all.vars(expr), valid_columns)
        if (length(invalid_cols) > 0) {
            if (verbose) {
                cat(sprintf("    [FAIL] Invalid columns found: %s\n", 
                    paste(invalid_cols, collapse = ", ")))
            }
            stop(sprintf(
                "Invalid columns in %s: %s",
                context_name, paste(invalid_cols, collapse = ", ")
            ))
        }
        if (verbose) cat("    [PASS] All columns valid\n")
    }
    
    tryCatch({
        # Check comparisons
        if (verbose) cat("Validating comparisons:\n")
        lapply(names(comparisons), function(comp_name) {
            check_expr_vars(comparisons[[comp_name]], 
                          sprintf("comparison '%s'", comp_name))
        })
        
        # Check control factors
        if (verbose) cat("Validating control factors:\n")
        lapply(names(control_factors), function(factor_name) {
            if (verbose) {
                cat(sprintf("  Checking control factor '%s'\n", factor_name))
            }
            invalid_cols <- setdiff(control_factors[[factor_name]], valid_columns)
            if (length(invalid_cols) > 0) {
                if (verbose) {
                    cat(sprintf("    [FAIL] Invalid columns: %s\n",
                        paste(invalid_cols, collapse = ", ")))
                }
                stop(sprintf(
                    "Invalid columns in control factor '%s': %s",
                    factor_name, paste(invalid_cols, collapse = ", ")
                ))
            }
            if (verbose) cat("    [PASS] All columns valid\n")
        })
        
        # Check conditions
        if (verbose) cat("Validating experimental conditions:\n")
        lapply(names(conditions), function(cond_name) {
            check_expr_vars(conditions[[cond_name]], 
                          sprintf("condition '%s'", cond_name))
        })
        
        if (verbose) cat("[PASS] All column references valid\n")
        return(TRUE)
    }, error = function(e) {
        if (verbose) cat(sprintf("[FAIL] Column reference validation: %s\n", e$message))
        if (stop_on_error) stop(e) else return(FALSE)
    })
}

validate_column_order <- function(
    categories,    # list: category definitions
    column_order,  # character: ordered column names
    stop_on_error = TRUE,
    verbose = FALSE
) {
    tryCatch({
        if (verbose) {
            cat("[VALIDATING] Checking column order...\n")
            cat("  Category columns:", paste(sort(names(categories)), collapse = ", "), "\n")
            cat("  Column order:", paste(sort(column_order), collapse = ", "), "\n")
        }
        
        if (!identical(sort(names(categories)), sort(column_order))) {
            if (verbose) cat("  [FAIL] Column order mismatch\n")
            stop("Column order must include all category columns")
        }
        
        if (verbose) cat("[PASS] Column order valid\n")
        return(TRUE)
    }, error = function(e) {
        if (verbose) cat(sprintf("[FAIL] Column order validation: %s\n", e$message))
        if (stop_on_error) stop(e) else return(FALSE)
    })
}
