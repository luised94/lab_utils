#!/usr/bin/env Rscript

library(curl)

# Function to download a file
download_file <- function(url, dest_path) {
    tryCatch({
        curl_download(url, dest_path, mode = "wb")
        cat(sprintf("Successfully downloaded: %s\n", dest_path))
    }, error = function(e) {
        stop(sprintf("Error downloading %s: %s", url, e$message))
    })
}

# Function to create directory if it doesn't exist
ensure_dir <- function(dir_path) {
    if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    cat(sprintf("Created directory: %s\n", dir_path))
    }
}

# Main function
main <- function() {
    base_dir <- file.path(Sys.getenv("HOME"), "data", "feature_files")
    ensure_dir(base_dir)
    
    cat("Defining the URLs and File paths...\n")
    # Define URLs and file paths
    data_sources <- list(
        hawkins = list(
            url = "https://ars.els-cdn.com/content/image/1-s2.0-S2211124713005834-mmc2.xlsx",
            path = file.path(base_dir, "hawkins-origins-timing.xlsx")
        ),

        eaton_peaks = list(
            url = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM424nnn/GSM424494/suppl/GSM424494_wt_G2_orc_chip_combined.bed.gz",
            path = file.path(base_dir, "eaton_peaks.bed.gz")
        ),

        sgd_features = list(
            url = "https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab",
            path = file.path(base_dir, "SGD_features.tab")
        ),

        sgd_gff = list(
            url = "https://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz",
            path = file.path(base_dir, "saccharomyces_cerevisiae.gff.gz")
        ),

        eaton_acs = list(
            url = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM424nnn/GSM424494/suppl/GSM424494_acs_locations.bed.gz",
            path = file.path(base_dir, "eaton_acs.bed.gz")
        )
    )
    
    
    # Download files
    cat("Starting data download...\n")
    for (source in names(data_sources)) {
        cat(sprintf("Downloading %s data...\n", source))
        download_file(data_sources[[source]]$url, data_sources[[source]]$path)
    }
    
    cat("All downloads complete.\n")
    
    # Optionally, add data validation here
    # For example, check file sizes or run basic parsing
}

# Run the main function
if (!interactive()) {
    main()
}
