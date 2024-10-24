#' Data Preparation Functions
prepare_visualization_data <- function(directory,
                                    chromosome = CONFIG$DEFAULTS$CHROMOSOME,
                                    genome_dir = CONFIG$DEFAULTS$GENOME_DIR,
                                    genome_pattern = CONFIG$DEFAULTS$GENOME_PATTERN) {
    log_info("Preparing visualization data")
    
    # Validate directory
    directory_path <- validate_input(directory)
    
    # Load sample table
    log_info("Loading sample table")
    sample_table <- load_sample_table(directory_path)
    
    # Load reference genome
    log_info("Loading reference genome")
    ref_genome <- load_reference_genome(
        genome_dir = genome_dir,
        genome_pattern = genome_pattern
    )
    
    # Create genome range
    log_info("Creating genome range")
    genome_range <- create_chromosome_GRange(ref_genome)
    
    list(
        directory = directory_path,
        samples = sample_table,
        genome = ref_genome,
        range = genome_range
    )
}

prepare_feature_track <- function(chromosome,
                                genome_range,
                                pattern = CONFIG$DEFAULTS$FEATURE_PATTERN) {
    log_info("Preparing feature track")
    
    feature_data <- load_feature_file_GRange(
        chromosome_to_plot = chromosome,
        feature_file_pattern = pattern,
        genomeRange_to_get = genome_range
    )
    
    AnnotationTrack(
        feature_data,
        name = paste("Origin Peaks", "Eaton 2010", sep = "")
    )
}
