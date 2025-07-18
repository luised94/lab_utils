################################################################################
# Flow cytometry Experiment Setup
# Author: Luis | Date: 2025-03-11 | Version: 1.0.0
################################################################################
# PURPOSE:
#   Creates standardized directory structure and metadata files for flow cytometry experiments
#
# USAGE:
# 1. Update experiment_id (format: YYMMDDBel, e.g., "241122Bel")
# 2. Run from cli.
#
# DEPENDENCIES: ~/lab_utils/core_scripts/configuration_experiment_bmc
#
# OUTPUTS:
# - Standard directory structure (peak/, fastq/, alignment/, bigwig/, plots/)
# - Sample metadata files (configuration_flow_cytometry.csv)
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

message("All function files loaded...")
################################################################################
# Handle script arguments
################################################################################
# Parse arguments and validate configurations
description <- "Setup flow cytometry experiments"
args <- parse_flow_cytometry_arguments(description = description)
args_info <- list(
    title = "Script Configuration",
    "script.name" = get_script_name(),
    "script.description" = description
)

print_debug_info(modifyList(args_info, args))
directory_path <- args$directory_path
experiment_id <- args$experiment_id
accept_configuration <- args$accept_configuration
experiment_dir <- args$experiment_dir

# Experiment ID Validation
stopifnot(
    "Only one experiment id required for this script" = length(experiment_id) == 1
)

message("Arguments parsed...")
################################################################################
# Load and Validate Experiment Configuration
################################################################################
# flow_cytometry_config.R is ignored in the git repository as this file is changed to add new experiments.
config_path <- "~/lab_utils/core_scripts/configuration_flow_cytometry.R"
# Define required dependencies
required_modules <- list(
    list(
        path = config_path,
        description = "Flow Cytometry Configuration",
        required = TRUE
    )
)

flow_cytometry_configuration_definition_path <- required_modules[[
    which(sapply(required_modules, function(x)
        x$description == "Flow Cytometry Configuration"
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
    ".total_modules" = length(required_modules),
    ".required_modules" = sum(sapply(required_modules, `[[`, "required"))
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

xit_file <- list.files(
    path = directory_path,
    pattern = paste0(experiment_id, "\\.xit$"),
    recursive = FALSE,
    include.dirs = FALSE
)

stopifnot(
    "Script experiment_id is not the same as CONFIG EXPERIMENT_ID" = experiment_id == EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID,
    "Only one xit file expected. Double check directory." = length(xit_file) == 1
)

message("Modules loaded...")
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

    print_debug_info(
        modifyList(
            list(
                title = "Final Configuration",
                "override.mode" = override_result$mode
                ),
            RUNTIME_CONFIG
        )
    )
    message("Override complete...")
}

handle_configuration_checkpoint(
    accept_configuration = accept_configuration,
    experiment_id = experiment_id
)

message("Configuration accepted...")
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
metadata <- metadata[
    do.call(
        order,
        metadata[EXPERIMENT_CONFIG$COLUMN_ORDER]
    ),
]

#--------------------
# Verify sample count
#--------------------
n_samples <- nrow(metadata)
expected <- EXPERIMENT_CONFIG$METADATA$EXPECTED_SAMPLES
if (n_samples != expected) {
    terminal_width <- get_width()
    # Set width to 200
    options(width = terminal_width)

    # Print diagnostic information
    cat("\nDiagnostic Information:\n")
    cat("----------------------\n")
    print(table(metadata[, category_to_show]))  # Show antibody distribution
    cat("\nFull sample breakdown:\n")
    print(summary(metadata))         # Show all category distributions
    cat("\n")

    # Control this in the configuration_flow_cytometry.R file.
    # Helps display metadata values to help narrow down where you need to add combinations to filter.
    if (show_all_metadata) {
        print(metadata)
    } else if (show_particular_metadata) {
        print(metadata[metadata[, category_to_show] == values_to_show, ])
    }

    stop(sprintf("Expected %d samples, got %d", expected, n_samples))
}

# Verification message
cat(sprintf("Total samples: %d\n", nrow(metadata)))
cat("\nDiagnostic Information:\n")
cat("----------------------\n")
print(table(metadata[, category_to_show]))
cat("----------------------\n")

# Generate sample names
metadata$full_name <- apply(metadata, 1, paste, collapse = "_")
metadata$short_name <- apply(metadata[, EXPERIMENT_CONFIG$COLUMN_ORDER], 1,
    function(x) paste0(substr(x, 1, 1), collapse = ""))

################################################################################
# File Output Generation
################################################################################
filenames <- c("flow_cytometry_config.R", "sample_grid.csv")
for (filename in filenames) {
    output_file_path <- file.path(directory_path, paste0(experiment_id, "_", filename))
    # Handle file writing with dry run checks
    if (RUNTIME_CONFIG$output_dry_run) {
        # Check if the file exists.
        if (file.exists(output_file_path)) {
            cat(sprintf("[DRY RUN] File exists: %s. Overwrite will occur if dry-run is disabled.\n", output_file_path))
            next
        }
        # Dry-run message
        cat(sprintf("[DRY RUN] Would write file to: %s\n", output_file_path))
        next
    }

    # Write the file. Overwrite by default.
    if(endsWith(filename, ".csv")) {
            safe_write_file(
                data = metadata,
                path = output_file_path,
                write_fn = write.csv,
                verbose = RUNTIME_CONFIG$debug_verbose,
                interactive = interactive(),
                row.names = FALSE
            )
        next
    }

    if(endsWith(filename, ".R")) {
        safe_write_file(
            data = flow_cytometry_configuration_definition_path,
            path = output_file_path,
            write_fn = file.copy,
            verbose = RUNTIME_CONFIG$debug_verbose,
            interactive = interactive(),
            overwrite = TRUE
        )
        next
    }

    warning(sprintf("Unsupported file extension for file: %s", filename))
}

print(metadata)
message("Setup for flow cytometry experiment complete...")
