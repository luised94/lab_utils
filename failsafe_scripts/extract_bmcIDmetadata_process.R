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


# Function to validate and enforce factor levels
enforce_factor_levels <- function(data_frame, categories) {
    print("Validating and enforcing factor levels")
    
    if (!all(names(categories) %in% colnames(data_frame))) {
        missing_cols <- setdiff(names(categories), colnames(data_frame))
        stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
    }
    
    for (col_name in names(categories)) {
        print(sprintf("Processing column: %s", col_name))
        
        # Check for invalid levels
        invalid_levels <- setdiff(data_frame[[col_name]], categories[[col_name]])
        if (length(invalid_levels) > 0) {
            stop(sprintf(
                "Invalid levels in %s: %s\nAllowed levels: %s",
                col_name,
                paste(invalid_levels, collapse = ", "),
                paste(categories[[col_name]], collapse = ", ")
            ))
        }
        
        # Convert to factor with predefined levels
        data_frame[[col_name]] <- factor(
            x = data_frame[[col_name]],
            levels = categories[[col_name]],
            ordered = TRUE
        )
        
        print(sprintf(
            "Enforced %d levels for %s",
            length(categories[[col_name]]),
            col_name
        ))
    }
    
    return(data_frame)
}

# Function to sort dataframe by multiple columns
sort_metadata_frame <- function(data_frame, column_order) {
    print("Sorting metadata frame")
    
    if (!all(column_order %in% colnames(data_frame))) {
        missing_cols <- setdiff(column_order, colnames(data_frame))
        stop(sprintf("Missing sort columns: %s", paste(missing_cols, collapse = ", ")))
    }
    
    sorted_frame <- data_frame[do.call(
        order,
        data_frame[column_order]
    ), ]
    
    print(sprintf("Sorted by columns: %s", paste(column_order, collapse = ", ")))
    
    return(sorted_frame)
}

# Modified process_metadata function
process_metadata <- function(directory_path, output_to_file = TRUE) {
    print("Starting metadata processing")
    
    # Source configuration
    source("~/lab_utils/scripts/bmc_config.R")
    if (!exists("EXPERIMENT_CONFIG")) {
        stop("Configuration loading failed: EXPERIMENT_CONFIG not found")
    }
    
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
    
    # Enforce factor levels
    metadata <- enforce_factor_levels(
        data_frame = metadata,
        categories = EXPERIMENT_CONFIG$CATEGORIES
    )
    
    # Sort dataframe
    metadata <- sort_metadata_frame(
        data_frame = metadata,
        column_order = EXPERIMENT_CONFIG$COLUMN_ORDER
    )
    
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
    metadata$sample_id <- exp_numbers
    
    # Save processed file
    if(output_to_file) {
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
    }
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

## Execute main if script is run directly
#if (identical(environment(), globalenv())) {
#    main()
#}
