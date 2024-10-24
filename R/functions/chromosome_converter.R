#' Chromosome Name Conversion Functions
convert_chromosome_names <- function(names,
                                   target_style,
                                   config = CONFIG$CHROMOSOMES) {
    log_info("Converting chromosome names to style:", target_style)
    
    vapply(
        names,
        function(chr) convert_single_chromosome(chr, target_style, config),
        character(1)
    )
}

convert_single_chromosome <- function(chr,
                                    target_style,
                                    config = CONFIG$CHROMOSOMES) {
    # Remove existing prefix if present
    chr <- gsub(CONFIG$NAMING$PATTERNS$PREFIX, "", chr, 
                ignore.case = TRUE)
    
    # Handle special chromosomes
    if (chr %in% config$SPECIAL) {
        return(paste0(config$PREFIX, chr))
    }
    
    # Convert based on current format
    if (grepl(CONFIG$NAMING$PATTERNS$ROMAN, chr)) {
        return(handle_roman_chromosome(chr, target_style))
    } else if (grepl(CONFIG$NAMING$PATTERNS$NUMERIC, chr)) {
        return(handle_numeric_chromosome(chr, target_style))
    }
    
    log_warning("Unable to convert chromosome:", chr)
    return(paste0(config$PREFIX, chr))
}

handle_roman_chromosome <- function(chr, target_style) {
    switch(target_style,
        "ENSEMBL" = as.character(arabic(chr)),
        "UCSC" = paste0(CONFIG$CHROMOSOMES$PREFIX, chr),
        "NCBI" = chr,
        stop("Unknown target style:", target_style)
    )
}

handle_numeric_chromosome <- function(chr, target_style) {
    switch(target_style,
        "ENSEMBL" = chr,
        "UCSC" = paste0(CONFIG$CHROMOSOMES$PREFIX, as.roman(as.integer(chr))),
        "NCBI" = as.roman(as.integer(chr)),
        stop("Unknown target style:", target_style)
    )
}
