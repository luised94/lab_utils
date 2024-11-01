#' BAM File Management Functions
find_sample_bam <- function(sample_info,
                           bam_dir,
                           reference_pattern = CONFIG$FILES$PATTERNS$REFERENCE) {
    log_info("Finding BAM file for:", sample_info$short_name)
    
    bam_files <- find_matching_bams(
        sample_info$sample_ID,
        bam_dir
    )
    
    if (length(bam_files) == 0) {
        log_warning("No BAM files found for:", sample_info$sample_ID)
        return(NULL)
    }
    
    # Filter for reference genome
    ref_files <- filter_reference_bams(
        bam_files,
        reference_pattern
    )
    
    if (length(ref_files) == 0) {
        log_warning("No reference BAM files found for:", sample_info$sample_ID)
        return(NULL)
    }
    
    ref_files[1]
}

find_matching_bams <- function(sample_id, bam_dir) {
    pattern <- sprintf(".*%s.*%s", sample_id, CONFIG$FILES$PATTERNS$BAM)
    list.files(bam_dir, pattern = pattern, full.names = TRUE)
}

filter_reference_bams <- function(bam_files,
                                reference_pattern) {
    bam_files[grepl(reference_pattern, bam_files)]
}

validate_bam_file <- function(file_path, sample_info) {
    if (!file.exists(file_path)) {
        log_warning(sprintf(
            "BAM file not found for sample %s (%s)",
            sample_info$short_name,
            sample_info$sample_ID
        ))
        return(FALSE)
    }
    return(TRUE)
}
