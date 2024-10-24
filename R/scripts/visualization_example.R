#!/usr/bin/env Rscript
# functions moved

source("../functions/track_generator.R")
source("../functions/range_handler.R")

#' Generate example visualization
create_example_visualization <- function(chromosome = CONFIG$GENOME$DEFAULT_CHROMOSOME,
                                      start = 1000,
                                      end = 5000) {
    log_info("Creating example visualization")
    
    # Create example data
    highlights <- create_example_highlights(chromosome, start, end)
    tracks <- create_example_tracks(chromosome, start, end)
    
    # Create highlight track
    highlight_track <- create_highlight_track(
        tracks,
        highlights,
        chromosome
    )
    
    # Generate plot
    plotTracks(
        highlight_track,
        from = start,
        to = end,
        chromosome = chromosome
    )
}

#' Example data generation functions
create_example_highlights <- function(chromosome, start, end) {
    GRanges(
        seqnames = chromosome,
        ranges = IRanges(
            start = c(start + 500, start + 2000),
            end = c(start + 1000, start + 2500)
        )
    )
}

create_example_tracks <- function(chromosome, start, end) {
    length <- end - start + 1
    
    list(
        create_genome_axis_track(chromosome),
        create_data_track(runif(length), chromosome, "Random"),
        create_data_track(rnorm(length), chromosome, "Normal"),
        create_annotation_track(chromosome, start, end)
    )
}
