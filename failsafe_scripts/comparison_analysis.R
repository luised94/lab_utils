# Function to load and validate metadata
load_metadata <- function(directory_path) {
    print("Loading processed metadata file")
    
    csv_path <- file.path(
        directory_path,
        "documentation",
        sprintf("%s_processed_grid.csv", basename(directory_path))
    )
    
    if (!file.exists(csv_path)) {
        stop(sprintf("Processed metadata file not found: %s", csv_path))
    }
    
    metadata <- read.csv(
        file = csv_path,
        stringsAsFactors = FALSE,
        check.names = TRUE
    )
    
    print(sprintf("Loaded metadata with %d rows", nrow(metadata)))
    return(metadata)
}

# Function to perform comparisons
analyze_comparisons <- function(metadata, comparisons) {
    print("Starting comparison analysis")
    
    results <- list()
    
    for (comp_name in names(comparisons)) {
        print(sprintf("\nAnalyzing comparison: %s", comp_name))
        
        # Evaluate the comparison expression
        subset_rows <- eval(comparisons[[comp_name]], metadata)
        subset_data <- metadata[subset_rows, ]
        
        # Store results
        results[[comp_name]] <- subset_data
        
        # Log results
        print(sprintf("Found %d matching samples:", nrow(subset_data)))
        if (nrow(subset_data) > 0) {
            print("Sample details:")
            for (i in 1:nrow(subset_data)) {
                print(sprintf(
                    "  Sample %d: %s_%s_%s_%s (Exp: %s)",
                    i,
                    subset_data$antibody[i],
                    subset_data$rescue_allele[i],
                    subset_data$auxin_treatment[i],
                    subset_data$time_after_release[i],
                    subset_data$experiment_number[i]
                ))
            }
        } else {
            print("  No matching samples found")
        }
    }
    
    return(invisible(results))
}

# Main execution
main <- function() {
    args <- commandArgs(trailingOnly = TRUE)
    
    if (length(args) != 1) {
        stop("Usage: Rscript comparison_analysis.R <experiment_directory>")
    }
    
    directory_path <- args[1]
    
    if (!dir.exists(directory_path)) {
        stop(sprintf("Directory not found: %s", directory_path))
    }
    
    print("Loading configuration")
    source("~/lab_utils/scripts/bmc_config.R")
    
    if (!exists("EXPERIMENT_CONFIG")) {
        stop("Configuration loading failed: EXPERIMENT_CONFIG not found")
    }
    
    print(sprintf("Analyzing experiment: %s", EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID))
    
    # Load metadata
    metadata <- load_metadata(directory_path)
    
    # Perform comparisons
    results <- analyze_comparisons(metadata, EXPERIMENT_CONFIG$COMPARISONS)
    
    print("\nAnalysis completed successfully")
}

# Execute main if script is run directly
if (identical(environment(), globalenv())) {
    main()
}
