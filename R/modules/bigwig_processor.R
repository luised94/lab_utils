#' BigWig File Processing Functions
find_bigwig_file <- function(directory,
                           sample_id,
                           pattern) {
    log_info("Finding bigwig file for:", sample_id)
    
    matches <- list.files(
        directory,
        pattern = sample_id,
        full.names = TRUE,
        recursive = TRUE
    )
    
    bigwig_files <- matches[grepl(pattern, matches)]
    
    if (length(bigwig_files) == 0) {
        log_warning("No bigwig file found for:", sample_id)
        return(NULL)
    }
    
    if (length(bigwig_files) > 1) {
        log_warning("Multiple bigwig files found, using first")
    }
    
    bigwig_files[1]
}

import_bigwig_data <- function(file_path, region) {
    log_info("Importing bigwig data from:", file_path)
    
    tryCatch({
        import(con = file_path, which = region)
    }, error = function(e) {
        log_error("Failed to import bigwig:", e$message)
        NULL
    })
}
