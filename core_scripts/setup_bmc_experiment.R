################################################################################
# BMC Experiment Setup Script
################################################################################
#
# PURPOSE:
#   Creates standardized directory structure and metadata files for BMC ChIP-seq
#   experiments, including sample tracking and submission documents.
#
# USAGE:
#   1. Use '/!!' in vim/neovim to jump to required updates
#   2. Update experiment_id (format: YYMMDD'Bel', e.g., "241122Bel")
#   3. Set RUNTIME_CONFIG options as needed
#   4. Source or run script
#
# !! ----> REQUIRED UPDATES:
# !! experiment_id <- "241122Bel"
# INPUTS:
#   - experiment_id: 9-character experiment identifier
#   - bmc_config.R: External configuration defining experimental design
#
# OUTPUTS:
#   1. Directory Structure:
#      ~/data/[experiment_id]/
#      +-- peak/
#      +-- fastq/
#      |   +-- raw/
#      |   +-- processed/
#      +-- alignment/
#      +-- bigwig/
#      +-- plots/
#      +-- documentation/
#
#   2. Files:
#      - [experiment_id]_sample_grid.csv: Complete experimental design
#      - [experiment_id]_bmc_table.tsv: BMC submission metadata
#      - [experiment_id]_bmc_config.R: Configuration snapshot
#
# CONTROLS:
#   RUNTIME_CONFIG$output_dry_run    = TRUE   # Preview without creating files
#   RUNTIME_CONFIG$debug_verbose    = TRUE   # Show detailed progress
#   RUNTIME_CONFIG$debug_interactive = TRUE  # Confirm before proceeding
# DEPENDENCIES:
#   - R base packages only
#   - ~/lab_utils/core_scripts/bmc_config.R
#
# COMMON ISSUES:
#   1. Wrong experiment ID format -> Check YYMMDD pattern
#   2. Unexpected sample count -> Review antibody distribution
#   3. File access denied -> Check ~/data permissions
#
# AUTHOR: Luis
# DATE: 2024-11-27
# VERSION: 2.0.0
#
################################################################################
################################################################################
source(file.path(Sys.getenv("HOME"), "lab_utils", "core_scripts", "functions_for_script_control.R"))

# Parse arguments and validate configurations
args <- parse_args(commandArgs(trailingOnly = TRUE))
experiment_id <- args[["experiment-id"]]
source(file.path("~/data", experiment_id, "documentation", 
                paste0(experiment_id, "_bmc_config.R")))
validate_configs(c("RUNTIME_CONFIG", "EXPERIMENT_CONFIG"))
################################################################################
# Experiment ID Validation
################################################################################
stopifnot(
    "Experiment ID must be a character string" = is.character(experiment_id),
    "Invalid experiment ID format. Expected: YYMMDD'Bel'" = grepl("^\\d{6}Bel$", experiment_id)
)

################################################################################
# Load and Validate Experiment Configuration
################################################################################
# Bootstrap phase
bootstrap_path <- normalizePath("~/lab_utils/core_scripts/functions_for_file_operations.R", 
                              mustWork = FALSE)
if (!file.exists(bootstrap_path)) {
    stop(sprintf("[FATAL] Bootstrap file not found: %s", bootstrap_path))
}
source(bootstrap_path)

# !! Update the path and the file accordingly.
config_path <- "~/bmc_config.R"
# Define required dependencies
required_modules <- list(
    list(
        path = config_path,
        description = "BMC Configuration",
        required = TRUE
    )
)

bmc_configuration_definition_path <- required_modules[[
    which(sapply(required_modules, function(x) 
        x$description == "BMC Configuration"
    ))
]]$path

# Validate module structure
stopifnot(
    "modules must have required fields" = all(sapply(required_modules, function(m) {
        all(c("path", "description", "required") %in% names(m))
    }))
)

# Load dependencies with status tracking
load_status <- lapply(required_modules, function(module) {
    if (RUNTIME_CONFIG$debug_verbose) {
        cat(sprintf("\n[LOADING] %s\n", module$description))
    }
    
    success <- safe_source(module$path, verbose = TRUE)
    
    if (!success && module$required) {
        stop(sprintf(
            "[FATAL] Failed to load required module: %s\n  Path: %s",
            module$description, module$path
        ))
    } else if (!success) {
        warning(sprintf(
            "[WARNING] Optional module not loaded: %s\n  Path: %s",
            module$description, module$path
        ))
    }
    
    return(list(
        module = module$description,
        path = module$path,
        loaded = success
    ))
})

# Display loading summary using ASCII
if (RUNTIME_CONFIG$debug_verbose) {
    cat("\n=== Module Loading Summary ===\n")
    invisible(lapply(load_status, function(status) {
        cat(sprintf(
            "[%s] %s\n    Path: %s\n",
            if(status$loaded) "+" else "-",
            status$module,
            status$path
        ))
    }))
}

stopifnot("Script experiment_id is not the same as CONFIG EXPERIMENT_ID" = experiment_id == EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID)

################################################################################
# Directory Setup and User Confirmation
################################################################################
base_dir <- file.path(Sys.getenv("HOME"), "data", experiment_id)
if (RUNTIME_CONFIG$debug_interactive) {
    cat(sprintf("\nExperiment ID: %s\n", experiment_id))
    cat("Base directory will be:", base_dir, "\n")

    user_response <- readline("Continue with this experiment? (y/n): ")
    if (tolower(user_response) != "y") {
        stop("Script terminated by user")
    }
    cat("Proceeding with directory creation...\n\n")
}

################################################################################
# Directory Structure Definition and Creation
################################################################################
# Define directory structure
data_directories <- c(
    "peak",
    "fastq/raw",
    "fastq/processed",
    "quality_control",
    "alignment",
    "bigwig",
    "plots/genome_tracks/overview",
    "plots/genome_tracks/experimental_comparisons",
    "documentation/dna_qc_traces",
    "documentation/config"
)

# Create directory structure
full_paths <- file.path(base_dir, data_directories)
invisible(lapply(full_paths, function(path) {
    if (RUNTIME_CONFIG$output_dry_run) {
        cat(sprintf("[DRY RUN] Would create directory: %s\n", path))
    } else {
        dir_created <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
        if (RUNTIME_CONFIG$debug_verbose) {
            status <- if (dir_created) "Created" else "Already exists"
            cat(sprintf("[%s] %s\n", status, path))
        }
    }
}))

# Report directory creation status
if (RUNTIME_CONFIG$debug_verbose) {
    mode <- if (RUNTIME_CONFIG$output_dry_run) "DRY RUN" else "LIVE RUN"
    cat(sprintf("\n[%s] Directory structure for experiment: %s\n", mode, experiment_id))
    cat(sprintf("[%s] Base directory: %s\n", mode, base_dir))
}

cat("Directories created successfully!\n")

################################################################################
# Sample Metadata Generation and Validation
################################################################################
# Generate experimental combinations
metadata <- do.call(expand.grid, EXPERIMENT_CONFIG$CATEGORIES)

# Filter invalid combinations
invalid_idx <- Reduce(
    `|`,
    lapply(EXPERIMENT_CONFIG$INVALID_COMBINATIONS, eval, envir = metadata)
)
metadata <- subset(metadata, !invalid_idx)

# Apply experimental conditions
#valid_idx <- Reduce(
#    `|`,
#    lapply(EXPERIMENT_CONFIG$EXPERIMENTAL_CONDITIONS, eval, envir = metadata)
#)
#metadata <- subset(metadata, valid_idx)

# Verify sample count
n_samples <- nrow(metadata)
expected <- EXPERIMENT_CONFIG$METADATA$EXPECTED_SAMPLES
if (n_samples != expected) {
    # Print diagnostic information
    cat("\nDiagnostic Information:\n")
    cat("----------------------\n")
    print(table(metadata$antibody))  # Show antibody distribution
    cat("\nFull sample breakdown:\n")
    print(summary(metadata))         # Show all category distributions
    cat("\n")

    stop(sprintf("Expected %d samples, got %d", expected, n_samples))
}

################################################################################
# Sample Classification
################################################################################
sample_classifications <- EXPERIMENT_CONFIG$SAMPLE_CLASSIFICATIONS

# First, create a matrix/data frame to store all classification results
classification_results <- matrix(FALSE, 
                               nrow = nrow(metadata), 
                               ncol = length(sample_classifications),
                               dimnames = list(NULL, names(sample_classifications)))

# Evaluate each classification condition
for (type in names(sample_classifications)) {
    classification_results[, type] <- eval(sample_classifications[[type]], 
                                         envir = metadata)
}

# Create the final classification vector
metadata$sample_type <- "treatment"  # Default classification
for (type in names(sample_classifications)) {
    # Find rows where this classification is TRUE
    matching_rows <- classification_results[, type]
    # Assign the type name (removing 'is_' prefix)
    metadata$sample_type[matching_rows] <- sub("^is_", "", type)
}

# Validation check
multiple_classifications <- rowSums(classification_results) > 1
if (any(multiple_classifications)) {
    cat("\nERROR: Multiple Classification Detected!\n")
    cat("----------------------------------------\n")
    
    # Show problematic samples with their classifications
    problem_samples <- metadata[multiple_classifications, ]
    cat("Samples with multiple classifications:\n\n")
    
    # Show which classifications were TRUE for each problematic sample
    for (i in which(multiple_classifications)) {
        cat(sprintf("\nSample %d:\n", i))
        cat("Sample details:\n")
        print(metadata[i, ])
        cat("\nMatching classifications:\n")
        matching_types <- names(classification_results[i,])[classification_results[i,]]
        print(matching_types)
        cat("----------------------------------------\n")
    }
    
    stop("Please fix multiple classifications in experiment configuration")
}

# Success diagnostic display
cat("\nSample Classification Summary:\n")
cat("============================\n")

# Overall counts
cat("\n1. Distribution of sample types:\n")
print(table(metadata$sample_type))

# Detailed breakdown by relevant factors
cat("\n2. Sample types by antibody:\n")
print(table(metadata$sample_type, metadata$antibody))

# Show a few samples from each classification
cat("\n3. Example samples from each classification:\n")
for (type in unique(metadata$sample_type)) {
    cat(sprintf("\n%s samples:\n", toupper(type)))
    print(metadata[metadata$sample_type == type, ][1:min(3, sum(metadata$sample_type == type)), ])
    cat("----------------------------------------\n")
}

# Verification message
cat("\nClassification Verification:\n")
cat(sprintf("- Total samples: %d\n", nrow(metadata)))
cat(sprintf("- Classified samples: %d\n", sum(table(metadata$sample_type))))
cat(sprintf("- Unclassified samples: %d\n", sum(is.na(metadata$sample_type))))

################################################################################
# Metadata Formatting and Organization
################################################################################
# Enforce factor levels from config
for (col_name in names(EXPERIMENT_CONFIG$CATEGORIES)) {
    if (col_name %in% colnames(metadata)) {
        metadata[[col_name]] <- factor(
            metadata[[col_name]],
            levels = EXPERIMENT_CONFIG$CATEGORIES[[col_name]],
            ordered = TRUE
        )
    }
}

# Sort metadata according to column order
metadata <- metadata[do.call(
    order,
    metadata[EXPERIMENT_CONFIG$COLUMN_ORDER]
), ]

# Generate sample names
metadata$full_name <- apply(metadata, 1, paste, collapse = "_")
metadata$short_name <- apply(metadata[, EXPERIMENT_CONFIG$COLUMN_ORDER], 1,
    function(x) paste0(substr(x, 1, 1), collapse = ""))

################################################################################
# BMC Metadata Generation
################################################################################
bmc_metadata <- data.frame(
    SampleName = metadata$full_name,
    Vol_uL = 10,
    Conc = 0,
    Type = "ChIP",
    Genome = "Saccharomyces cerevisiae",
    Notes = ifelse(
        metadata$antibody == "Input",
        "Run on fragment analyzer.",
        "Run on femto pulse."
    ),
    Pool = "A",
    stringsAsFactors = FALSE
)

################################################################################
# File Output Generation
################################################################################
# Define output file paths
sample_grid_path <- file.path(base_dir, "documentation",
                             paste0(experiment_id,"_", "sample_grid.csv"))

bmc_table_path <- file.path(base_dir, "documentation",
                           paste0(experiment_id,"_", "bmc_table.tsv"))

bmc_experiment_config_path <- file.path(base_dir, "documentation",
                           paste0(experiment_id,"_", "bmc_config.R"))

# Handle file writing with dry run checks
if (RUNTIME_CONFIG$output_dry_run) {
    cat(sprintf("[DRY RUN] Would write sample grid to: %s\n", sample_grid_path))
    cat(sprintf("[DRY RUN] Would write BMC table to: %s\n", bmc_table_path))
    cat(sprintf("[DRY RUN] Would write BMC config script: %s\n", bmc_experiment_config_path))
} else {
    # Write sample grid file
    safe_write_file(
        data = metadata,
        path = sample_grid_path,
        write_fn = write.csv,
        verbose = RUNTIME_CONFIG$debug_verbose,
        row.names = FALSE
    )

    # Write BMC table file
    safe_write_file(
        data = bmc_metadata,
        path = bmc_table_path,
        write_fn = write.table,
        verbose = RUNTIME_CONFIG$debug_verbose,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    # Copy BMC Experiment Config
    safe_write_file(
        data = bmc_configuration_definition_path,
        path = bmc_experiment_config_path,
        write_fn = file.copy,
        verbose = RUNTIME_CONFIG$debug_verbose,
        overwrite = TRUE
    )
}
