#' Data Download Functions
download_feature_data <- function(source_config,
                                base_dir,
                                timestamp) {
    log_info("Downloading feature data")
    
    results <- list()
    
    for (source in names(source_config)) {
        result <- download_source_data(
            source,
            source_config[[source]],
            base_dir,
            timestamp
        )
        
        results[[source]] <- result
    }
    
    results
}

download_source_data <- function(source_name,
                               config,
                               base_dir,
                               timestamp) {
    log_info("Downloading:", source_name)
    
    output_file <- generate_output_path(
        base_dir,
        config$FILE,
        timestamp
    )
    
    # Download file
    download_file(config$URL, output_file)
    
    # Process if compressed
    if (is_compressed(output_file)) {
        output_file <- decompress_file(output_file)
    }
    
    # Validate and process
    process_downloaded_file(
        output_file,
        source_name,
        config$TYPE
    )
}

download_file <- function(url, dest_path) {
    log_info("Downloading:", basename(dest_path))
    
    tryCatch({
        curl::curl_download(url, dest_path, mode = "wb")
        log_info("Download complete")
        TRUE
    }, error = function(e) {
        log_error("Download failed:", e$message)
        FALSE
    })
}
