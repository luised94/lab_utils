#' Mapping Statistics Functions
calculate_mapping_stats <- function(flagstat_files,
                                  sample_info,
                                  genome_names) {
    log_info("Calculating mapping statistics")
    
    # Initialize results matrix
    results <- initialize_results_matrix(genome_names)
    
    # Process each sample
    sample_ids <- get_sample_ids(sample_info)
    
    for (sample_id in sample_ids) {
        stats <- process_sample_stats(
            sample_id,
            flagstat_files,
            genome_names
        )
        
        results <- rbind(results, stats)
    }
    
    results
}

process_sample_stats <- function(sample_id,
                               flagstat_files,
                               genome_names) {
    log_info("Processing sample:", sample_id)
    
    # Find relevant flagstat files
    sample_files <- find_sample_flagstats(
        sample_id,
        flagstat_files
    )
    
    # Calculate percentages
    stats <- sapply(sample_files, function(file) {
        calculate_mapping_percentage(file)
    })
    
    # Create result row
    result <- as.data.frame(t(stats))
    colnames(result) <- genome_names
    
    result
}

calculate_mapping_percentage <- function(flagstat_file) {
    log_info("Calculating mapping percentage for:", 
             basename(flagstat_file))
    
    # Read flagstat data
    data <- read.table(flagstat_file, sep = "\t")
    
    # Find relevant rows
    patterns <- CONFIG$BAM_QC$METRICS
    is_relevant <- Reduce("|", lapply(patterns, grepl, data[,3]))
    subset_data <- data[is_relevant, ]
    
    # Calculate percentage
    mapped <- as.numeric(subset_data[2,1])
    total <- as.numeric(subset_data[1,1])
    
    (mapped / total) * 100
}
