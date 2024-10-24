#' Main Plot Management Function
plot_all_sample_tracks <- function(sample_table,
                                 directory_path,
                                 chromosome_to_plot = 10,
                                 genomeRange_to_get,
                                 control_track,
                                 annotation_track,
                                 highlight_gr) {
    log_info("Starting sample track plotting")
    
    # Setup
    setup_result <- setup_plotting_environment(
        directory_path,
        chromosome_to_plot
    )
    
    # Process samples
    for (sample_index in seq_len(nrow(sample_table))) {
        if (sample_index == 1) {  # Only process first sample for now
            process_sample(
                sample_table[sample_index, ],
                sample_table,
                setup_result,
                control_track,
                annotation_track,
                highlight_gr
            )
        }
    }
    
    log_info("Plotting completed")
}

#' Helper Functions
setup_plotting_environment <- function(directory_path,
                                     chromosome) {
    log_info("Setting up plotting environment")
    
    list(
        plot_dir = file.path(directory_path, CONFIG$PATHS$SUBDIRS$PLOTS),
        bigwig_dir = file.path(directory_path, CONFIG$PATHS$SUBDIRS$BIGWIG),
        chromosome = paste0("chr", as.roman(chromosome)),
        title = sprintf("Complete View of Chrom %s", chromosome)
    )
}

process_sample <- function(sample_row,
                         sample_table,
                         setup,
                         control_track,
                         annotation_track,
                         highlight_gr) {
    log_info("Processing sample:", sample_row$short_name)
    
    # Get sample data
    sample_data <- get_sample_data(
        sample_row,
        setup$bigwig_dir,
        genomeRange_to_get
    )
    
    if (is.null(sample_data)) {
        log_warning("No data for sample:", sample_row$short_name)
        return(NULL)
    }
    
    # Get control data
    control_data <- get_control_data(
        sample_row,
        sample_table,
        setup$bigwig_dir,
        genomeRange_to_get
    )
    
    # Create tracks
    tracks <- create_track_set(
        sample_data,
        control_data,
        annotation_track,
        highlight_gr,
        setup$chromosome
    )
    
    # Generate plot
    output_file <- generate_output_filename(
        setup$plot_dir,
        sample_row$short_name,
        setup$chromosome
    )
    
    generate_sample_plot(
        tracks,
        setup$chromosome,
        setup$title,
        output_file
    )
}
