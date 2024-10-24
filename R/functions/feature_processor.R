#' Feature File Processing Functions
load_feature_file_GRange <- function(chromosome_to_plot = 10,
                                   feature_file_pattern = CONFIG$FEATURE_TYPES$PEAKS,
                                   genomeRange_to_get) {
    log_info("Loading feature file:", feature_file_pattern)
    
    # Validate directory
    feature_dir <- CONFIG$PATHS$FEATURE_DIR
    if (!dir.exists(feature_dir)) {
        log_error("Feature directory not found:", feature_dir)
        stop("Missing feature directory")
    }
    
    # Find feature file
    feature_files <- list.files(feature_dir,
                              pattern = feature_file_pattern,
                              full.names = TRUE,
                              recursive = TRUE)
    
    if (length(feature_files) != 1) {
        log_error("Invalid number of feature files:", length(feature_files))
        stop("Feature file error")
    }
    
    # Load and process features
    feature_grange <- tryCatch({
        import.bed(feature_files[1])
    }, error = function(e) {
        log_error("Failed to import feature file:", e$message)
        stop("Import error")
    })
    
    # Process chromosome styles
    feature_style <- determine_chr_style(seqlevels(feature_grange))
    genome_style <- determine_chr_style(seqlevels(genomeRange_to_get))
    
    log_info("Feature style:", feature_style)
    log_info("Genome style:", genome_style)
    
    # Normalize chromosome styles
    feature_grange_subset <- normalize_feature_range(
        feature_grange,
        genomeRange_to_get,
        feature_style,
        genome_style
    )
    
    return(feature_grange_subset)
}

normalize_feature_range <- function(feature_grange,
                                  genome_range,
                                  feature_style,
                                  genome_style) {
    if (feature_style == genome_style) {
        log_info("Styles match, using direct subsetting")
        return(subsetByOverlaps(feature_grange, genome_range))
    }
    
    log_info("Adjusting styles for compatibility")
    
    # Adjust genome range to match feature style
    adjusted_genome_range <- genome_range
    new_seqlevels <- normalize_chr_names(seqlevels(genome_range),
                                       feature_style)
    seqlevels(adjusted_genome_range) <- new_seqlevels
    
    # Subset features
    feature_subset <- subsetByOverlaps(feature_grange,
                                     adjusted_genome_range)
    
    # Convert back to genome style
    final_seqlevels <- normalize_chr_names(seqlevels(feature_subset),
                                         genome_style)
    seqlevels(feature_subset) <- final_seqlevels
    
    return(feature_subset)
}
#' Feature File Processing Functions
process_feature_files <- function(input_dir = CONFIG$FEATURES$PATHS$BASE_DIR) {
    log_info("Processing feature files")
    
    # Validate directory
    validate_directory(input_dir)
    
    # Get file list
    files <- get_feature_files(input_dir)
    
    if (length(files) == 0) {
        log_warning("No files to process")
        return(NULL)
    }
    
    # Process each file
    results <- process_files(files)
    
    # Verify outputs
    verify_outputs(input_dir)
    
    results
}

process_files <- function(files) {
    log_info("Processing", length(files), "files")
    
    lapply(files, function(file) {
        tryCatch({
            process_single_file(file)
        }, error = function(e) {
            log_error("Failed to process:", basename(file))
            log_error("Error:", e$message)
            NULL
        })
    })
}

process_single_file <- function(file_path) {
    log_info("Processing file:", basename(file_path))
    
    # Read data
    data <- read_feature_file(file_path)
    
    # Process data
    processed <- process_feature_data(data, basename(file_path))
    
    # Convert to GRanges
    granges <- convert_to_granges(processed, basename(file_path))
    
    # Output results
    output_results(granges, file_path)
    
    granges
}
