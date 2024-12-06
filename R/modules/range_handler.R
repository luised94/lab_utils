#' Genome Range Management Functions
create_chromosome_range <- function(genome_data,
                                  config = CONFIG$GENOME) {
    log_info("Creating chromosome range")
    
    validate_genome_data(genome_data)
    
    GRanges(
        seqnames = genome_data$chrom,
        ranges = IRanges(
            start = 1,
            end = genome_data$basePairSize
        ),
        strand = config$DEFAULT_STRAND
    )
}

load_control_range <- function(control_dir,
                             identifier,
                             chromosome,
                             genome_range) {
    log_info("Loading control range data")
    
    bigwig_file <- find_control_bigwig(
        control_dir,
        identifier
    )
    
    if (is.null(bigwig_file)) {
        log_error("Control bigwig file not found")
        return(NULL)
    }
    
    import_control_data(
        bigwig_file,
        chromosome,
        genome_range
    )
}

import_control_data <- function(file_path,
                              chromosome,
                              genome_range) {
    log_info("Importing control data")
    
    control_style <- determine_chr_style(
        seqlevels(import(file_path))
    )
    
    chromosome_name <- normalize_chr_names(
        chromosome,
        control_style
    )
    
    subset_range <- genome_range[
        seqnames(genome_range) == chromosome_name
    ]
    
    import(file_path, which = subset_range)
}
