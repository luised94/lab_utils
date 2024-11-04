# Function to extract experiment numbers from fastq files
extract_experiment_numbers <- function(directory_path) {
    print(sprintf("Scanning directory: %s", directory_path))
    
    fastq_files <- list.files(
        path = file.path(directory_path, "fastq"),
        pattern = "consolidated_.*_sequence\\.fastq$",
        full.names = FALSE
    )
    
    if (length(fastq_files) == 0) {
        stop("No fastq files found in specified directory")
    }
    
    print(sprintf("Found %d fastq files", length(fastq_files)))
    
    # Extract numbers using base R
    exp_numbers <- gsub(
        pattern = "consolidated_([0-9]{5,6})_sequence\\.fastq",
        replacement = "\\1",
        x = fastq_files
    )
    
    return(exp_numbers)
}

# Function to process metadata and add experiment numbers
process_metadata <- function(directory_path) {
    print("Starting metadata processing")
    
    # Construct file path
    csv_path <- file.path(
        directory_path,
        "documentation",
        sprintf("%s_sample_grid.csv", basename(directory_path))
    )
    
    # Validate file existence
    if (!file.exists(csv_path)) {
        stop(sprintf("Metadata file not found: %s", csv_path))
    }
    
    # Read CSV
    metadata <- read.csv(
        file = csv_path,
        stringsAsFactors = FALSE,
        check.names = TRUE
    )
    
    print(sprintf("Read metadata file with %d rows", nrow(metadata)))
    
    # Get experiment numbers
    exp_numbers <- extract_experiment_numbers(directory_path)
    
    # Validate length match
    if (length(exp_numbers) != nrow(metadata)) {
        stop(sprintf(
            "Mismatch between number of fastq files (%d) and metadata rows (%d)",
            length(exp_numbers),
            nrow(metadata)
        ))
    }
    
    # Add experiment numbers
    metadata$experiment_number <- exp_numbers
    
    # Save processed file
    output_path <- file.path(
        directory_path,
        "documentation",
        sprintf("%s_processed_grid.csv", basename(directory_path))
    )
    
    write.csv(
        x = metadata,
        file = output_path,
        row.names = FALSE,
        quote = TRUE
    )
    
    print(sprintf("Saved processed metadata to: %s", output_path))
    
    return(metadata)
}

# Main execution
main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    
    if (length(args) != 1) {
        stop("Usage: Rscript experiment_metadata_process.R <experiment_directory>")
    }
    
    directory_path <- args[1]
    
    if (!dir.exists(directory_path)) {
        stop(sprintf("Directory not found: %s", directory_path))
    }
    
    processed_data <- process_metadata(directory_path)
    print("Processing completed successfully")
}

# Execute main if script is run directly
if (identical(environment(), globalenv())) {
    main()
}
