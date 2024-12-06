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
#' Plot Generation Functions
generate_sample_plot <- function(tracks,
                               chromosome,
                               title,
                               output_file = NULL,
                               config = CONFIG$VISUALIZATION) {
    log_info("Generating plot:", title)
    
    if (!is.null(output_file)) {
        svg(output_file)
        on.exit(dev.off())
    }
    
    plotTracks(
        tracks,
        main = title,
        chromosome = chromosome,
        ylim = c(config$LIMITS$Y_MIN, config$LIMITS$Y_MAX)
    )
}

generate_output_filename <- function(base_dir,
                                   sample_name,
                                   chromosome,
                                   config = CONFIG$VISUALIZATION) {
    log_info("Generating output filename")
    
    date_str <- format(Sys.time(), config$OUTPUT$DATE_FORMAT)
    
    file.path(
        base_dir,
        sprintf(
            "%s_%s_%s_WithInputAndHighlights.%s",
            date_str,
            chromosome,
            sample_name,
            config$OUTPUT$FORMAT
        )
    )
}
