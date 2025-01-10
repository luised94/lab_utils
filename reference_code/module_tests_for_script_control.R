#!/usr/bin/env Rscript
source(file.path(Sys.getenv("HOME"), "lab_utils", "core_scripts", "functions_for_script_control.R"))
################################################################################
# Load Required Libraries
################################################################################
required_packages <- c("rtracklayer", "GenomicRanges", "Gviz")
check_required_packages(required_packages, verbose = TRUE)

################################################################################
# Handle script arguments
################################################################################
# Parse arguments and validate configurations
description <- "Script for testing all modules."
args <- parse_common_arguments(description = description)
experiment_id <- args$experiment_id
accept_configuration <- args$accept_configuration
experiment_dir <- args$experiment_dir
if (args$is_template) {
    config_path <- file.path(args$experiment_dir, "template_bmc_config.R")
    metadata_path <- file.path(args$experiment_dir, "template_sample_grid.csv")
} else {
    config_path <- file.path(args$experiment_dir, "documentation", paste0(args$experiment_id, "_bmc_config.R"))
    metadata_path <- file.path(args$experiment_dir, "documentation", paste0(args$experiment_id, "_sample_grid.csv"))
}
#config_path <- file.path(experiment_dir, "documentation", paste0(experiment_id, "_bmc_config.R"))
#metadata_path <- file.path(experiment_dir, "documentation", paste0(experiment_id, "_sample_grid.csv"))


################################################################################
# Load and Validate Experiment Configuration and Dependencies
################################################################################
# Bootstrap phase
bootstrap_path <- normalizePath("~/lab_utils/core_scripts/functions_for_file_operations.R",
                              mustWork = FALSE)
if (!file.exists(bootstrap_path)) {
    stop(sprintf("[FATAL] Bootstrap file not found: %s", bootstrap_path))
}
source(bootstrap_path)
safe_source(config_path)

# Define required dependencies
required_modules <- list(
    list(
        path = "~/lab_utils/core_scripts/functions_for_logging.R",
        description = "Logging functions",
        required = TRUE
    ),
    list(
        path = "~/lab_utils/core_scripts/functions_for_metadata_processing.R",
        description = "Process metadata grid for downstream analysis.",
        required = TRUE
    )
)

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

required_configs <- c("RUNTIME_CONFIG", "EXPERIMENT_CONFIG", "GENOME_TRACK_CONFIG")
validate_configs(required_configs)
invisible(lapply(required_configs, function(config) {
    print_config_settings(get(config), title = config)
}))

handle_configuration_checkpoint(
    accept_configuration = accept_configuration,
    experiment_id = experiment_id
)

# Setup directories
# !! Add directories here.
required_directories <- c("fastq", "documentation")
dirs <- setup_experiment_dirs(experiment_dir = experiment_dir,
    output_dir_name = "plots",
    required_input_dirs = required_directories
)
dirs$output_dir <- file.path(dirs$plots, "genome_tracks", "overview")
dir.create(dirs$output_dir, recursive = TRUE, showWarnings = FALSE)



# Find fastq files and extract sample IDs
fastq_files <- list.files(
    path = dirs$fastq,
    pattern = GENOME_TRACK_CONFIG$file_pattern,
    full.names = FALSE
)

if (length(fastq_files) == 0) {
    stop("No fastq files found in specified directory")
}

# Extract sample IDs from fastq filenames
sample_ids <- gsub(
    pattern = GENOME_TRACK_CONFIG$file_sample_id,
    replacement = "\\1",
    x = fastq_files
)

metadata <- load_and_process_experiment_metadata(
    metadata_path = metadata_path,
    categories = EXPERIMENT_CONFIG$CATEGORIES,
    column_order = EXPERIMENT_CONFIG$COLUMN_ORDER,
    sample_ids = sample_ids,
    stop_on_missing_columns = TRUE
)

if (RUNTIME_CONFIG$debug_verbose) {
    message("\nExperiment Setup Summary:")
    message(paste(rep("-", 50), collapse = ""))
    #message(sprintf("Input files found  : %d", length(fastq_files)))
    message(sprintf("Metadata rows      : %d", nrow(metadata)))
    message(sprintf("Output directory   : %s", dirs$output_dir))
    message(sprintf("Factor columns     : %s", 
                   paste(names(EXPERIMENT_CONFIG$CATEGORIES), collapse = ", ")))
    message(paste(rep("-", 50), collapse = ""))
    message(paste(rep("-", 50), collapse = ""))
    head(metadata)
    message(paste(rep("-", 50), collapse = ""))
}

