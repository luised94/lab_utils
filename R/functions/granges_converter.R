#' GRanges Conversion Functions
convert_granges_style <- function(gr,
                                target_style,
                                verbose = FALSE) {
    log_info("Converting GRanges object to style:", target_style)
    
    validate_granges(gr)
    
    # Get current seqnames
    current_seqnames <- as.character(seqnames(gr))
    if (verbose) {
        log_info("Original seqnames:",
                paste(unique(current_seqnames), collapse = ", "))
    }
    
    # Convert names
    new_seqnames <- convert_chromosome_names(
        current_seqnames,
        target_style
    )
    
    if (verbose) {
        log_info("Converted seqnames:",
                paste(unique(new_seqnames), collapse = ", "))
        log_info("Changes made:",
                sum(new_seqnames != current_seqnames))
    }
    
    # Update GRanges object
    update_granges_seqnames(gr, new_seqnames)
}

validate_granges <- function(gr) {
    if (!is(gr, "GRanges")) {
        log_error("Invalid input type")
        stop("Input must be a GenomicRanges object")
    }
}

update_granges_seqnames <- function(gr, new_names) {
    seqlevels(gr) <- unique(new_names)
    seqnames(gr) <- new_names
    gr
}
