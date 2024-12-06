#STATUS: REMOVE.
#!/usr/bin/env Rscript
library(readr)
library(dplyr)
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

# Adjust SGD_features column names
adjust_sgd_features <- function(df) {
    # Names taken from SGD_features.README
    new_names <- c(
        "Primary_SGDID", "Feature_type", "Feature_qualifier", "Feature_name",
        "Standard_gene_name", "Alias", "Parent_feature_name", "Secondary_SGDID",
        "Chromosome", "Start_coordinate", "Stop_coordinate", "Strand",
        "Genetic_position", "Coordinate_version", "Sequence_version",
        "Description"
    )

    if(ncol(df) != length(new_names)) {
        cat(sprintf("Number of columns in df: %s\nLength of name array: %s\n", ncol(df), length(new_names)))
        stop("Provided dataframe or table does not have same length as column names\n")
    } else {
        names(df) <- new_names
    }
    return(df)
}

# Adjust eaton acs bed file coordinates for proper loading as GRange object
adjust_eaton_acs <- function(file_path) {
    df <- read.table(file_path, sep = "\t", header = FALSE, col.names=c("chrom", "start", "end", "name", "score", "strand"))

    # Identify problematic rows
    df <- df %>%
        mutate(width = end - start,
               is_valid = width >= -1)

    problematic <- df %>% filter(!is_valid)
    cat(sprintf("Number of problematic rows %s:\n", nrow(problematic)))
    
    if (nrow(problematic) > 0) {
    # Fix the problematic rows by swapping start and end
        df_fixed <- df %>%
                    mutate(temp = ifelse(start > end, start, end),
                           start = ifelse(start > end, end, start),
                           end = temp) %>%
                    select(-temp, -width, -is_valid)
        
        fixed_file_path <- sub("\\.bed$", "_fixed.bed", file_path)
        write.table(df_fixed, fixed_file_path, sep ="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
        cat(sprintf("Fixed file written to:%s\n", fixed_file_path))
    } else {
        cat("No problematic rows found. File is valid.\n")
    }
}

# Function that orchestrates file readers based on extension
read_file <- function(file_path) {
    cat(sprintf("Reading file: %s\n", basename(file_path)))
    file_extension <- tools::file_ext(file_path)
    file_readers <- list(
        xls = readxl::read_excel,
        xlsx = readxl::read_excel,
        bed = rtracklayer::import.bed,
        gff3 = function(file_path) rtracklayer::import(file_path, format = 'gff3'),
        gff = function(file_path) rtracklayer::import(file_path, format = 'gff'),
        tab = function(file_path) read_tsv(file_path, col_names = FALSE),
        rds = readRDS
    )

    if (!(file_extension %in% names(file_readers))) {
        stop(sprintf("Unsupported file type: %s", file_extension))
    }
    data <- file_readers[[file_extension]](file_path)
    return(data)
}
# Main function
main <- function() {
    timestamp <- format(Sys.time(), "%y%m%d")
    base_dir <- file.path(Sys.getenv("HOME"), "data", "feature_files")
    ensure_dir(base_dir)
    
    cat("Defining the URLs and File paths...\n")
    # Define URLs and file paths
    data_sources <- list(
        hawkins = list(
            url = "https://ars.els-cdn.com/content/image/1-s2.0-S2211124713005834-mmc2.xlsx",
            path = "hawkins-origins-timing.xlsx"
        ),

        eaton_peaks = list(
            url = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM424nnn/GSM424494/suppl/GSM424494_wt_G2_orc_chip_combined.bed.gz",
            path = "eaton_peaks.bed.gz"
        ),

        sgd_features = list(
            url = "https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab",
            path = "SGD_features.tab"
        ),

        sgd_gff = list(
            url = "https://downloads.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz",
            path = "saccharomyces_cerevisiae.gff.gz"
        ),

        eaton_acs = list(
            url = "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM424nnn/GSM424494/suppl/GSM424494_acs_locations.bed.gz",
            path = "eaton_acs.bed.gz"
        )
    )
    
    
    # Download files
    cat("Starting data download...\n")
    for (source in names(data_sources)) {
        cat(sprintf("Downloading %s data...\n", source))
        output_file <- file.path(base_dir , paste0(timestamp, "_", data_sources[[source]]$path))
        download_file(data_sources[[source]]$url, output_file)
        if (grepl("\\.gz$", output_file)) {
            cat(sprintf("Gunzip %s file\n", output_file))
            R.utils::gunzip(output_file, remove = TRUE)
        }
        
        cat(sprintf("Download %s complete...\n", source))
        cat(sprintf("Verifying %s data...\n", source))
        output_file <- gsub("\\.gz$", "", output_file)

        if (source == "sgd_features") {
            sgd_df <- read_file(output_file)
            sgd_df <- adjust_sgd_features(sgd_df)
            cat(sprintf("\nHead of %s:\n", source))
            print(head(sgd_df))
        } else if (source == "eaton_acs") {
            fixed_file_path <- sub("\\.bed$", "_fixed.bed", output_file)
            adjust_eaton_acs(output_file)
            df <- read_file(fixed_file_path)
            cat(sprintf("\nHead of %s:\n", source))
            print(head(df))
        } else {
            df <- read_file(output_file)
            cat(sprintf("\nHead of %s:\n", source))
            print(head(df))
        }

    }
    
    cat("Verified all files\n")
    cat("rm 240830_eaton_acs.bed if you dont need the file.")
    cat("Done.")

    
    # Optionally, add data validation here
    # For example, check file sizes or run basic parsing
}

# Run the main function
if (!interactive()) {
    main()
}
