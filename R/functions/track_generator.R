#' Track Generation Functions
create_genome_axis_track <- function(chromosome, name = NULL) {
    log_info("Creating genome axis track")
    
    GenomeAxisTrack(
        name = name %||% sprintf("Chr %s Axis", chromosome)
    )
}

create_data_track <- function(data,
                            chromosome,
                            name,
                            config = CONFIG$VISUALIZATION) {
    log_info("Creating data track:", name)
    
    if (!is(data, "GRanges")) {
        data <- convert_to_granges(data, chromosome)
    }
    
    DataTrack(
        data,
        type = config$TRACKS$TYPES$DATA,
        name = name,
        col = config$TRACKS$COLORS$PRIMARY,
        chromosome = chromosome
    )
}

create_highlight_track <- function(track_list,
                                 highlights,
                                 chromosome) {
    log_info("Creating highlight track")
    
    if (!is(highlights, "GRanges")) {
        stop("Highlights must be a GRanges object")
    }
    
    HighlightTrack(
        trackList = track_list,
        start = start(highlights),
        end = end(highlights),
        chromosome = as.character(seqnames(highlights))
    )
}
