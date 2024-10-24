#' Track Assembly Functions
create_track_set <- function(sample_data,
                           control_data,
                           annotation_data,
                           highlight_data,
                           chromosome,
                           config = CONFIG$VISUALIZATION) {
    log_info("Creating track set")
    
    # Create base tracks
    tracks <- list(
        create_genome_axis_track(chromosome),
        create_control_track(control_data, chromosome),
        create_sample_track(sample_data, chromosome),
        create_annotation_track(annotation_data, chromosome)
    )
    
    # Add highlights if present
    if (!is.null(highlight_data)) {
        tracks <- create_highlight_track(
            tracks,
            highlight_data,
            chromosome
        )
    }
    
    return(tracks)
}

create_sample_track <- function(data,
                              chromosome,
                              name,
                              config = CONFIG$VISUALIZATION) {
    log_info("Creating sample track:", name)
    
    DataTrack(
        data,
        type = config$TRACK_TYPES$LINE,
        name = name,
        col = config$COLORS$SAMPLE,
        chromosome = chromosome
    )
}

create_control_track <- function(data,
                               chromosome,
                               name,
                               config = CONFIG$VISUALIZATION) {
    log_info("Creating control track:", name)
    
    DataTrack(
        data,
        type = config$TRACK_TYPES$LINE,
        name = name,
        col = config$COLORS$CONTROL,
        chromosome = chromosome
    )
}
