#' Plot Generation Functions
generate_track_plot <- function(tracks,
                              title,
                              chromosome,
                              output_file = NULL) {
    log_info("Generating track plot:", title)
    
    if (!is.null(output_file)) {
        svg(output_file)
        on.exit(dev.off())
    }
    
    plotTracks(
        tracks,
        main = title,
        chromosome = chromosome
    )
}

create_output_filename <- function(base_dir,
                                 chromosome,
                                 pattern,
                                 comparison,
                                 time_id = NULL) {
    if (is.null(time_id)) {
        time_id <- format(Sys.time(), CONFIG$PLOT$DATE_FORMAT)
    }
    
    pattern_clean <- gsub("_|\\.bw", "", pattern)
    comparison_clean <- gsub("_", "", comparison)
    
    file.path(
        base_dir,
        sprintf(
            "%s_%s_%s_%s_%s.svg",
            time_id,
            chromosome,
            pattern_clean,
            comparison_clean,
            format(Sys.time(), "%H%M%S")
        )
    )
}

