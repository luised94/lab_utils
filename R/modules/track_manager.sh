#' Track Creation and Management Functions
create_genome_track <- function(chromosome) {
    log_info("Creating genome axis track")
    
    GenomeAxisTrack(
        name = sprintf("Chr %s Axis", chromosome)
    )
}

create_data_track <- function(bigwig_data,
                            name,
                            chromosome,
                            config = CONFIG$PLOT) {
    log_info("Creating data track:", name)
    
    DataTrack(
        bigwig_data,
        type = config$TRACK_TYPE,
        name = name,
        col = config$COLORS$DEFAULT,
        chromosome = chromosome
    )
}

create_control_track <- function(bigwig_data,
                               chromosome,
                               config = CONFIG$PLOT) {
    log_info("Creating control track")
    
    DataTrack(
        bigwig_data,
        type = config$TRACK_TYPE,
        name = "Input",
        col = config$COLORS$CONTROL,
        chromosome = chromosome
    )
}
