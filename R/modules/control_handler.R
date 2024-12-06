#' Control Sample Management Functions
find_control_sample <- function(sample_row,
                              sample_table,
                              factors,
                              bigwig_dir,
                              pattern) {
    log_info("Finding control sample")
    
    control_index <- determine_matching_control(
        sample_row,
        sample_table,
        factors
    )
    
    if (length(control_index) == 0) {
        log_warning("No matching control found, using fallback")
        return(find_fallback_control(sample_table, bigwig_dir, pattern))
    }
    
    control_data <- get_control_data(
        sample_table[control_index, ],
        bigwig_dir,
        pattern
    )
    
    if (is.null(control_data)) {
        log_warning("Control data not found, using fallback")
        return(find_fallback_control(sample_table, bigwig_dir, pattern))
    }
    
    control_data
}

get_control_data <- function(control_sample,
                           bigwig_dir,
                           pattern) {
    log_info("Getting control data for:", control_sample$sample_ID)
    
    bigwig_file <- find_bigwig_file(
        bigwig_dir,
        control_sample$sample_ID,
        pattern
    )
    
    if (is.null(bigwig_file)) {
        return(NULL)
    }
    
    list(
        id = control_sample$sample_ID,
        name = control_sample$short_name,
        file = bigwig_file
    )
}
