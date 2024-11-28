validate_category_values <- function(
    categories,     # list: category name -> character vector mapping
    stop_on_error = TRUE  # logical: whether to stop or return validation status
) {
    tryCatch({
        sapply(names(categories), function(category_name) {
            values <- categories[[category_name]]
            stopifnot(
                sprintf("Category '%s' values must be character vectors", category_name) = 
                    is.character(values),
                sprintf("Category '%s' values must be unique", category_name) = 
                    !any(duplicated(values))
            )
        })
        return(TRUE)
    }, error = function(e) {
        if (stop_on_error) stop(e) else return(FALSE)
    })
}

validate_column_references <- function(
    categories,        # list: valid column names
    comparisons,       # list: named expressions
    control_factors,   # list: factor definitions
    conditions,        # list: experimental conditions
    stop_on_error = TRUE
) {
    valid_columns <- names(categories)
    
    # Helper function for expression validation
    check_expr_vars <- function(expr, context_name) {
        invalid_cols <- setdiff(all.vars(expr), valid_columns)
        if (length(invalid_cols) > 0) {
            stop(sprintf(
                "Invalid columns in %s: %s",
                context_name, paste(invalid_cols, collapse = ", ")
            ))
        }
    }
    
    tryCatch({
        # Check comparisons
        lapply(names(comparisons), function(comp_name) {
            check_expr_vars(comparisons[[comp_name]], 
                          sprintf("comparison '%s'", comp_name))
        })
        
        # Check control factors
        lapply(names(control_factors), function(factor_name) {
            invalid_cols <- setdiff(control_factors[[factor_name]], valid_columns)
            if (length(invalid_cols) > 0) {
                stop(sprintf(
                    "Invalid columns in control factor '%s': %s",
                    factor_name, paste(invalid_cols, collapse = ", ")
                ))
            }
        })
        
        # Check conditions
        lapply(names(conditions), function(cond_name) {
            check_expr_vars(conditions[[cond_name]], 
                          sprintf("condition '%s'", cond_name))
        })
        
        return(TRUE)
    }, error = function(e) {
        if (stop_on_error) stop(e) else return(FALSE)
    })
}

validate_column_order <- function(
    categories,    # list: category definitions
    column_order,  # character: ordered column names
    stop_on_error = TRUE
) {
    tryCatch({
        stopifnot(
            "Column order must include all category columns" =
                identical(sort(names(categories)), sort(column_order))
        )
        return(TRUE)
    }, error = function(e) {
        if (stop_on_error) stop(e) else return(FALSE)
    })
}
