#' Enhanced control sample management
find_valid_control <- function(sample_row,
                             sample_table,
                             bam_dir,
                             factors) {
    log_info("Finding valid control sample")
    
    # Find matching control
    control_index <- determine_matching_control(
        sample_row,
        sample_table,
        factors
    )
    
    # Validate and possibly adjust control
    control_index <- validate_control_index(
        control_index,
        sample_table,
        bam_dir
    )
    
    if (control_index > 0) {
        return(get_control_bam(
            sample_table[control_index, ],
            bam_dir
        ))
    }
    
    log_warning("No valid control found")
    return(NULL)
}

validate_control_index <- function(index,
                                 sample_table,
                                 bam_dir) {
    if (length(index) == 0) {
        log_warning("No matching control found, using fallback")
        return(find_fallback_control(sample_table, bam_dir))
    }
    
    if (length(index) > CONFIG$CONTROL$MAX_CONTROLS) {
        log_warning("Multiple controls found, using first")
        index <- index[1]
    }
    
    # Validate BAM file exists
    control_bam <- find_matching_bam(
        sample_table$sample_ID[index],
        bam_dir
    )
    
    if (!is.null(control_bam)) {
        return(index)
    }
    
    log_warning("Control BAM not found, searching for alternative")
    return(find_fallback_control(sample_table, bam_dir))
}

find_fallback_control <- function(sample_table, bam_dir) {
    log_info("Searching for fallback control")
    
    input_samples <- sample_table[
        sample_table[[CONFIG$CONTROL$ANTIBODY_COLUMN]] == 
        CONFIG$CONTROL$INPUT_VALUE,
    ]
    
    for (i in seq_len(nrow(input_samples))) {
        if (!is.null(find_matching_bam(
            input_samples$sample_ID[i],
            bam_dir
        ))) {
            return(which(sample_table$sample_ID == 
                        input_samples$sample_ID[i]))
        }
    }
    
    return(CONFIG$CONTROL$DEFAULT_INDEX)
}
