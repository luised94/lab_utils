#' Data Conversion Functions
convert_to_granges <- function(data,
                             file_name,
                             config = CONFIG$FEATURES) {
    log_info("Converting data to GRanges")
    
    if (is(data, "GRanges")) {
        return(data)
    }
    
    converter <- determine_converter(file_name)
    converter(data, file_name)
}

determine_converter <- function(file_name) {
    for (type in names(CONFIG$FEATURES$FILE_TYPES)) {
        if (grepl(CONFIG$FEATURES$FILE_TYPES[[type]], file_name)) {
            return(get(paste0("convert_", tolower(type))))
        }
    }
    
    function(data, file_name) {
        log_warning("No specific converter for:", file_name)
        as(data, "GRanges")
    }
}

convert_nucleosome <- function(data, file_name) {
    log_info("Converting nucleosome data")
    
    config <- CONFIG$FEATURES$COLUMNS$NUCLEOSOME
    
    # Create metadata
    metadata <- data[, !(colnames(data) %in% config$EXCLUDE)]
    colnames(metadata) <- clean_column_names(colnames(metadata))
    
    GRanges(
        seqnames = data$`Nucleosome ID`,
        ranges = IRanges(
            start = data[[config$POSITION]],
            end = data[[config$POSITION]]
        ),
        strand = "*",
        chromosome = data$Chromosome,
        metadata
    )
}

convert_timing <- function(data, file_name) {
    log_info("Converting timing data")
    
    config <- CONFIG$FEATURES$COLUMNS$TIMING
    window <- config$WINDOW
    
    # Create origin names
    names <- paste0(data$Chromosome, "_", data$Position)
    
    # Create metadata
    metadata <- data[, !(colnames(data) %in% config$EXCLUDE)]
    colnames(metadata) <- clean_column_names(colnames(metadata))
    
    GRanges(
        seqnames = names,
        ranges = IRanges(
            start = data$Position - window,
            end = data$Position + window
        ),
        strand = "*",
        chromosome = data$Chromosome,
        metadata
    )
}
