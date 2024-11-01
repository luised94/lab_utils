#' BAM QC Analysis Functions
analyze_bam_qc <- function(directory,
                          config = CONFIG$BAM_QC) {
    log_info("Analyzing BAM QC for:", directory)
    
    # Get flagstat files
    flagstat_files <- find_flagstat_files(directory)
    
    if (length(flagstat_files) == 0) {
        log_warning("No flagstat files found")
        return(NULL)
    }
    
    # Get genome names
    genome_names <- extract_genome_names(flagstat_files)
    
    # Get sample information
    sample_info <- load_sample_info(directory)
    
    # Calculate mapping percentages
    mapping_stats <- calculate_mapping_stats(
        flagstat_files,
        sample_info,
        genome_names
    )
    
    mapping_stats
}

find_flagstat_files <- function(directory) {
    pattern <- CONFIG$BAM_QC$PATTERNS$FLAGSTAT
    qc_dir <- file.path(directory, CONFIG$PATHS$QC_DIR)
    
    list.files(
        qc_dir,
        pattern = pattern,
        recursive = FALSE,
        full.names = TRUE
    )
}

extract_genome_names <- function(files) {
    unique(sapply(files, function(path) {
        parts <- strsplit(basename(path), "_")[[1]]
        parts[length(parts) - 1]
    }))
}

load_sample_info <- function(directory) {
    pattern <- CONFIG$BAM_QC$PATTERNS$SAMPLE_INFO
    doc_dir <- file.path(directory, CONFIG$PATHS$DOC_DIR)
    
    info_file <- list.files(
        doc_dir,
        pattern = pattern,
        recursive = FALSE,
        full.names = TRUE
    )[1]
    
    if (is.na(info_file)) {
        log_error("Sample info file not found")
        return(NULL)
    }
    
    read.table(info_file, sep = ",", header = TRUE)
}
