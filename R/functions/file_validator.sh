#' Merge with existing validation functions
validate_bam_file <- function(file_path, sample_info) {
    log_info("Validating BAM file:", file_path)
    
    if (!file.exists(file_path)) {
        log_warning(sprintf("File not found for sample %s (%s)",
                          sample_info$short_name,
                          sample_info$sample_ID))
        return(FALSE)
    }
    
    return(TRUE)
}

find_matching_bam <- function(sample_id,
                            bam_dir,
                            reference_pattern = CONFIG$PATTERNS$REFERENCE_GENOME) {
    log_info("Finding BAM file for:", sample_id)
    
    pattern <- paste0(".*", sample_id, ".*", CONFIG$PATTERNS$BAM_SUFFIX)
    all_files <- list.files(bam_dir, pattern = pattern, full.names = TRUE)
    
    ref_files <- all_files[grepl(reference_pattern, all_files)]
    
    if (length(ref_files) == 0) {
        log_warning("No matching BAM file found")
        return(NULL)
    }
    
    ref_files[1]
}
