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

# End message
message("Script completed succesfully...")
