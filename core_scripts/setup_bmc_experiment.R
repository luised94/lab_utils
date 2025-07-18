#!/usr/bin/env Rscript
################################################################################
# BMC ChIP-seq Experiment Setup
# Author: Luis | Date: 2024-11-27 | Version: 2.0.0
################################################################################
# PURPOSE:
#   Creates standardized directory structure and metadata files for BMC ChIP-seq experiments
#
# USAGE:
#   Update configuration_experiment_bmc file.
#   $./setup_bmc_experiment.R --experiment-id=<experiment-id> --accept-configuration
#
# DEPENDENCIES: ~/lab_utils/core_scripts/configuration_experiment_bmc
#
# OUTPUTS:
# - Standard directory structure (peak/, fastq/, alignment/, bigwig/, plots/)
# - Sample metadata files (_sample_grid.csv, _bmc_table.tsv)
################################################################################
# Bootstrap phase
function_filenames <- c("logging", "script_control", "file_operations")
for (function_filename in function_filenames) {
    function_filepath <- sprintf("~/lab_utils/core_scripts/functions_for_%s.R", function_filename)
    normalized_path <- normalizePath(function_filepath)
    if (!file.exists(normalized_path)) {
        stop(sprintf("[FATAL] File with functions not found: %s", normalized_path))
    }
    source(normalized_path)
}

################################################################################
# Handle script arguments
################################################################################
# Parse arguments and validate configurations
description <- "Setup experiment directory."
args <- parse_common_arguments(description = description)
experiment_id <- args$experiment_id
accept_configuration <- args$accept_configuration
experiment_dir <- args$experiment_dir

args_info <- list(
    title = "Script Configuration",
    "script.name" = get_script_name(),
    "script.description" = description
)
print_debug_info(modifyList(args_info, args))
################################################################################
# Experiment ID Validation
################################################################################
stopifnot(
    "Only one experiment id required for this script" = length(experiment_id) == 1
)

################################################################################
# Load and Validate Experiment Configuration
################################################################################
# !! Update the path and the file accordingly.
# configuration_experiment_bmc is ignored in the git repository as this file is changed to add new experiments.
config_path <- "~/lab_utils/core_scripts/configuration_experiment_bmc"
# Define required dependencies
required_modules <- list(
    list(
        path = config_path,
        description = "BMC Configuration",
        required = TRUE
    )
)

bmc_configuration_data_path <- required_modules[[
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

load_status <- lapply(required_modules, function(module) {
    success <- safe_source(module$path, verbose = TRUE)
    if (!success && module$required) {
        stop(sprintf("Failed to load required module: %s\n  Path: %s",
            module$description, module$path))
    } else if (!success) {
        warning(sprintf("Optional module not loaded: %s\n  Path: %s",
            module$description, module$path))
    }
    list(
        module = module$description,
        path = module$path,
        loaded = success,
        required = module$required
    )
})

# Create debug info structure
module_info <- list(
    title = "Module Loading Status",
    "total_modules" = length(required_modules),
    "required_modules" = sum(sapply(required_modules, `[[`, "required"))
)

for (status in load_status) {
    module_key <- paste0(
        if(status$required) "required." else "optional.",
        gsub(" ", "_", tolower(status$module))
    )
    module_info[[module_key]] <- sprintf(
        "%s (%s)",
        status$module,  # Now showing description
        if(status$loaded) sprintf("loaded from %s", status$path) else "failed"
    )
}
# Display using print_debug_info
print_debug_info(module_info)

required_configs <- c("EXPERIMENT_CONFIG", "RUNTIME_CONFIG")
validate_configs(required_configs)
invisible(lapply(required_configs, function(config) {
    print_config_settings(get(config), title = config)
}))
stopifnot("Script experiment_id is not the same as CONFIG EXPERIMENT_ID" = experiment_id == EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID)

################################################################################
# Directory Setup and User Confirmation
################################################################################
# Handle configuration override (independent)
if (!is.null(args$override)) {
    structured_log_info("Starting override")
    override_result <- apply_runtime_override(
        config = RUNTIME_CONFIG,
        preset_name = args$override,
        preset_list = OVERRIDE_PRESETS
    )
    RUNTIME_CONFIG <- override_result$modified

 print_debug_info(modifyList(
     list(
         title = "Final Configuration",
         "override.mode" = override_result$mode
     ),
     RUNTIME_CONFIG  # Flat list of current settings
 ))
}

handle_configuration_checkpoint(
    accept_configuration = accept_configuration,
    experiment_id = experiment_id
)

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
    "coverage",
    "plots/genome_tracks/overview",
    "plots/genome_tracks/experimental_comparisons",
    "documentation/dna_qc_traces",
    "documentation/config"
)

# Create directory structure
full_paths <- file.path(experiment_dir, data_directories)
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
    cat(sprintf("[%s] Base directory: %s\n", mode, experiment_dir))
}

cat("Directories created successfully!\n")

################################################################################
# Sample Metadata Generation and Validation
################################################################################
# Generate experimental combinations
metadata <- do.call(expand.grid, EXPERIMENT_CONFIG$CATEGORIES)

# Filter invalid combinations if set.
# Run without combination to see outputs.
if (length(EXPERIMENT_CONFIG$INVALID_COMBINATIONS) > 0) {
    invalid_idx <- Reduce(
        `|`,
        lapply(EXPERIMENT_CONFIG$INVALID_COMBINATIONS, eval, envir = metadata)
    )
    metadata <- subset(metadata, !invalid_idx)
}

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
    terminal_width <- get_width()
    # Set width to 200
    options(width = terminal_width)

    # Print diagnostic information
    cat("\nDiagnostic Information:\n")
    cat("----------------------\n")
    print(table(metadata$antibody))  # Show antibody distribution
    cat("\nFull sample breakdown:\n")
    print(summary(metadata))         # Show all category distributions
    cat("\n")

    # Control this in the configuration_experiment_bmc file.
    # Helps display metadata values to help narrow down where you need to add combinations to filter.
    if (show_all_metadata) {
        print(metadata)
    } else if (show_particular_metadata) {
        print(metadata[metadata[, category_to_show] == values_to_show, ])
    }

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
# Default classification
metadata$sample_type <- "treatment"
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
    Type = ifelse(
        metadata$antibody == "Input",
        "Input",
        "ChIP"
    ),
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
filenames <- c("sample_grid.csv", "bmc_table.tsv", "configuration_experiment_bmc")
# Loop through each filename to handle path assignment and file writing
for (filename in filenames) {
    # Construct the output file path
    output_file_path <- file.path(experiment_dir, "documentation", paste0(experiment_id, "_", filename))
    # Handle file writing with dry run checks
    if (RUNTIME_CONFIG$output_dry_run) {
        # Dry-run message
        if (file.exists(output_file_path)) {
            cat(sprintf("[DRY RUN] File exists: %s. Overwrite will occur if dry-run is disabled.\n", output_file_path))
        } else {
            cat(sprintf("[DRY RUN] Would write file to: %s\n", output_file_path))
        }
    } else {
        # Determine the appropriate write function based on the file extension
        if (endsWith(filename, ".csv")) {
            safe_write_file(
                data = metadata,
                path = output_file_path,
                write_fn = write.csv,
                verbose = RUNTIME_CONFIG$debug_verbose,
                interactive = interactive(),  # Use interactive() to detect if running interactively
                row.names = FALSE
            )
        } else if (endsWith(filename, ".tsv")) {
            safe_write_file(
                data = bmc_metadata,
                path = output_file_path,
                write_fn = write.table,
                verbose = RUNTIME_CONFIG$debug_verbose,
                interactive = interactive(),
                sep = "\t",
                row.names = FALSE,
                quote = FALSE
            )
        } else if (endsWith(filename, ".R")) {
            safe_write_file(
                data = bmc_configuration_data_path,
                path = output_file_path,
                write_fn = file.copy,
                verbose = RUNTIME_CONFIG$debug_verbose,
                interactive = interactive(),
                overwrite = TRUE
            )
        } else {
            warning(sprintf("Unsupported file extension for file: %s", filename))
        }
    }
}
print(metadata)
