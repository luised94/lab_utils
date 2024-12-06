#' Sample and Control Matching Functions
find_matching_samples <- function(sample_table,
                                directory_path,
                                config = CONFIG) {
    log_info("Finding matching samples and controls")
    
    # Setup directories
    directories <- setup_directories(directory_path)
    
    # Get matching factors
    factors <- get_factors_to_match(sample_table)
    
    # Process each sample
    results <- process_all_samples(
        sample_table,
        directories,
        factors
    )
    
    return(results)
}

setup_directories <- function(base_path) {
    list(
        bam = file.path(base_path, CONFIG$PATHS$SUBDIRS$ALIGNMENT),
        bigwig = file.path(base_path, CONFIG$PATHS$SUBDIRS$BIGWIG)
    )
}

process_all_samples <- function(sample_table,
                              directories,
                              factors) {
    log_info("Processing all samples")
    
    results <- list()
    
    for (i in seq_len(nrow(sample_table))) {
        result <- process_single_sample(
            sample_table[i, ],
            sample_table,
            directories,
            factors
        )
        
        if (!is.null(result)) {
            results[[i]] <- result
        }
    }
    
    return(results)
}

process_single_sample <- function(sample_row,
                                sample_table,
                                directories,
                                factors) {
    log_info("Processing sample:", sample_row$short_name)
    
    # Find control
    control_info <- find_control_sample(
        sample_row,
        sample_table,
        directories$bam,
        factors
    )
    
    if (is.null(control_info)) {
        log_warning("No valid control found for:", sample_row$short_name)
        return(NULL)
    }
    
    # Find sample BAM
    sample_bam <- find_sample_bam(
        sample_row,
        directories$bam
    )
    
    if (is.null(sample_bam)) {
        log_warning("No BAM file found for:", sample_row$short_name)
        return(NULL)
    }
    
    list(
        sample = sample_bam,
        control = control_info$file
    )
}
