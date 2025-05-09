#!/usr/bin/env Rscript
###############################################################################
# Plot bigwig files generated from bam from same sample but differentially processed
################################################################################
# PURPOSE: Plot bigwig files to see if blacklist filtering helps samples with many reads
#          at the ends or in blacklisted regions
# USAGE: ./plot_genome_tracks_of_different_processed_bam_files.R --experiment-id=<experiment-id> <options>
# DEPENDENCIES: GenomicRanges, rtracklayer
# OUTPUT: svg plots with comparisons for blacklist and non-blacklisted files.
# AUTHOR: LEMR
# DATE: 2025-05-09
################################################################################
current_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
# Bootstrap phase
# Also loads OVERRIDE_PRESETS
FUNCTION_FILENAMES <- c("logging", "script_control", "file_operations")
for (function_filename in FUNCTION_FILENAMES) {
    function_filepath <- sprintf("~/lab_utils/core_scripts/functions_for_%s.R", function_filename)
    normalized_path <- normalizePath(function_filepath)
    if (!file.exists(normalized_path)) {
        stop(sprintf("[FATAL] File with functions not found: %s", normalized_path))
    }
    source(normalized_path)
}
message("Bootstrap phase completed...")

################################################################################
# Handle script arguments
################################################################################
# Parse arguments and validate configurations
description <- "Plot genome tracks in batches"
args <- parse_common_arguments(description = description)
experiment_id <- args$experiment_id
accept_configuration <- args$accept_configuration
experiment_dir <- args$experiment_dir
is_template <- args$is_template

file_directory <- if (is_template) args$experiment_dir else file.path(args$experiment_dir, "documentation")
file_identifier <- if (is_template) "template" else args$experiment_id

config_path <- file.path(file_directory, paste0(file_identifier, "_bmc_config.R"))
metadata_path <- file.path(file_directory, paste0(file_identifier, "_sample_grid.csv"))
message("Arguments parsed...")

args_info <- list(
    title = "Script Configuration",
    "script.name" = get_script_name(),
    "script.description" = description
)
print_debug_info(modifyList(args_info, args))

################################################################################
# Load Required Libraries
################################################################################
required_packages <- c("rtracklayer", "GenomicRanges", "Gviz")
check_required_packages(
  packages = required_packages, verbose = TRUE,
  skip_validation = args$skip_validation
)

message("Packages confirmed...")
################################################################################
# Load and Validate Experiment Configuration and Dependencies
################################################################################
# TODO: Consolidate this into the previous for loop? Or use second for loop with safe_source //
# Would need to remove the print_debug_info call. I probably need to adjust that function anyways. //
# Define required dependencies
required_modules <- list(
    list(
        path = "~/lab_utils/core_scripts/functions_for_metadata_processing.R",
        description = "Process metadata grid for downstream analysis.",
        required = TRUE
    ),
    list(
        path = "~/lab_utils/core_scripts/functions_for_genome_tracks.R",
        description = "Functions to load genome track objects for plotting",
        required = TRUE
    ),
    list(
        path = config_path,
        description = "Functions to load genome track objects for plotting",
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
# Process module loading
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

# Add status for each module
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

required_configs <- c("EXPERIMENT_CONFIG", "GENOME_TRACK_CONFIG", "RUNTIME_CONFIG")
validate_configs(required_configs)
invisible(lapply(required_configs, function(config) {
    print_config_settings(get(config), title = config)
}))

if(!is.null(args$output_format)) {
    RUNTIME_CONFIG$output_format <- args$output_format

    message(sprintf("Setting output format: %s", RUNTIME_CONFIG$output_format))
}

# Handle configuration override (independent)
if (!is.null(args$override)) {
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

# End message -------
message("Script completed succesfully...")
