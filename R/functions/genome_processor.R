#' Genome Loading Functions
load_reference_genome <- function(genome_dir = CONFIG$PATHS$GENOME_DIR, 
                                genome_pattern = CONFIG$PATTERNS$GENOME) {
    log_info("Loading reference genome")
    
    directory_path <- file.path(CONFIG$PATHS$BASE_DIR, genome_dir)
    
    if (!dir.exists(directory_path)) {
        log_error("Genome directory not found:", directory_path)
        stop("Missing genome directory")
    }
    
    genome_files <- list.files(directory_path, pattern = genome_pattern, 
                             full.names = TRUE, recursive = TRUE)
    
    if (length(genome_files) != 1) {
        log_error("Invalid number of genome files found:", length(genome_files))
        stop("Genome file error")
    }
    
    genome <- readFasta(genome_files[1])
    
    genome_df <- data.frame(
        chrom = names(as(genome, "DNAStringSet")),
        basePairSize = width(genome)
    ) %>% 
        filter(chrom != "chrM")
    
    return(genome_df)
}

#' Chromosome Name Processing Functions
normalize_chr_names <- function(chr_names, target_style) {
    log_info("Normalizing chromosome names")
    
    if (!target_style %in% CONFIG$CHROMOSOME_MAPPING$STYLES) {
        log_error("Invalid chromosome style:", target_style)
        stop("Invalid style")
    }
    
    chr_names <- gsub(paste0("^", CONFIG$CHROMOSOME_MAPPING$PREFIX), "", 
                     chr_names)
    
    normalized <- switch(
        target_style,
        "UCSC" = paste0(CONFIG$CHROMOSOME_MAPPING$PREFIX, chr_names),
        "Roman" = paste0(CONFIG$CHROMOSOME_MAPPING$PREFIX, 
                        mapvalues(chr_names, 
                                from = as.character(1:16),
                                to = CONFIG$CHROMOSOME_MAPPING$ROMAN)),
        "Numeric" = mapvalues(chr_names,
                            from = CONFIG$CHROMOSOME_MAPPING$ROMAN,
                            to = as.character(1:16))
    )
    
    return(unname(normalized))
}
