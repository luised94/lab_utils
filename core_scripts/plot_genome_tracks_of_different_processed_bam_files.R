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

# Checkpoint handler ---------------------------
# See the settings before running. Override with --accept-configuration option.
handle_configuration_checkpoint(
  accept_configuration = accept_configuration,
  experiment_id = experiment_id
)

################################################################################
# Setup directories, variables and metadata
################################################################################
# !! Add directories here.
required_directories <- c("fastq", "documentation", "coverage")
# Creates directories and stores them in dirs variable.
dirs <- setup_experiment_dirs(
  experiment_dir = experiment_dir,
  output_dir_name = "plots",
  required_input_dirs = required_directories
)

# Control output directory here.
dirs$output_dir <- file.path(dirs$plots, "genome_tracks", "processing_intermediates")
dir.create(dirs$output_dir, recursive = TRUE, showWarnings = FALSE)

# Find fastq files and extract sample IDs
fastq_files <- list.files(
    path = dirs$fastq,
    pattern = GENOME_TRACK_CONFIG$file_pattern,
    full.names = TRUE
)

# Find bigwig files
bigwig_pattern <- sprintf(
    "processed_.*_sequence_to_S288C_%s\\.bw$",
    EXPERIMENT_CONFIG$NORMALIZATION$active
)

bigwig_files <- list.files(
    dirs$coverage,
    pattern = bigwig_pattern,
    full.names = TRUE
)

bigwig_basenames <- basename(bigwig_files)
fastq_basenames <- basename(fastq_files)

normalization_method <- sub(
    ".*_([^_]+)\\.bw$",
    "\\1",
    basename(bigwig_files[1])
)

if (length(bigwig_files) == 0) {
    stop("No bigwig files found in specified directory")
}

if (length(fastq_files) == 0) {
    warning("No fastq files found in specified directory")
    message("Attempting to use bigwig files to find sample_ids")
    # Extract sample IDs from fastq filenames
    sample_ids <- gsub(
        pattern = GENOME_TRACK_CONFIG$file_sample_id_from_bigwig,
        replacement = "\\1",
        x = bigwig_basenames
    )
    stopifnot(
        "Length of samples_ids is not lower than length of bigwig files." =
         all(nchar(sample_ids) < nchar(bigwig_basenames))
    )
} else {
    # Extract sample IDs from fastq filenames
    sample_ids <- gsub(
        pattern = GENOME_TRACK_CONFIG$file_sample_id,
        replacement = "\\1",
        x = fastq_basenames
    )
    stopifnot(
        "Length of samples_ids is not lower than length of fastq files." =
        all(nchar(sample_ids) < nchar(fastq_basenames))
    )
}

if (any(nchar(sample_ids) == 0)) {
  empty_ids <- bigwig_basenames[nchar(sample_ids) == 0]
  message(
    "Some sample_ids are empty.\n",
    "Number of files with missing ids:", length(empty_ids), "\n",
    "Affected files: ", paste(empty_ids, collapse = ", "), "\n"
    )
  stop("Some sample ids are empty.\nCheck your file name extraction pattern.\nEach sample id must be non-empty.")

}

if (RUNTIME_CONFIG$debug_verbose) {
  message("\nFile Processing:")
  message("  Directories:")
  message(sprintf("    FASTQ: %s", dirs$fastq))
  message(sprintf("    Coverage: %s", dirs$coverage))
  message("\n  Pattern Matching:")
  message(sprintf("    FASTQ pattern: %s", GENOME_TRACK_CONFIG$file_pattern))
  message(sprintf("    Bigwig pattern: %s", bigwig_pattern))
  message("\n  File Discovery:")
  message(sprintf("    FASTQ files found: %d", length(fastq_files)))
  if (length(fastq_files) > 0) {
    message("    First few FASTQ files:")
    invisible(lapply(head(fastq_files, 3), function(f) message("      ", f)))
  }
  message(sprintf("\n    Bigwig files found: %d", length(bigwig_files)))
  if (length(bigwig_files) > 0) {
    message("    First few Bigwig files:")
    invisible(lapply(head(bigwig_files, 3), function(f) message("      ", basename(f))))
    message(sprintf("\n    Normalization: %s", normalization_method))
  }
  message("\n  Sample IDs:")
  message(sprintf("    Total found: %d", length(sample_ids)))
  message("    First few IDs:")
  invisible(lapply(head(sample_ids, 3), function(id) message("      ", id)))
}

# Load csv, turn into factors, ensure proper order, add sample id column.
# INQ: Should be "cached"?
metadata <- load_and_process_experiment_metadata(
  metadata_path = metadata_path,
  categories = EXPERIMENT_CONFIG$CATEGORIES,
  column_order = EXPERIMENT_CONFIG$COLUMN_ORDER,
  sample_ids = sample_ids,
  stop_on_missing_columns = TRUE
)

# Create color scheme
color_scheme <- create_color_scheme(
    config = list(
        placeholder = GENOME_TRACK_CONFIG$color_placeholder,
        input = GENOME_TRACK_CONFIG$color_input
    ),
    categories = list(
        antibody = unique(metadata$antibody),
        rescue_allele = unique(metadata$rescue_allele)
    )
)

# Shorten sample ids to minimal distinguishing set
short_sample_ids <- create_minimal_identifiers(
  metadata$sample_id,
  verbose = RUNTIME_CONFIG$debug_verbose
)

# Create mapping between full and short IDs
sample_id_mapping <- setNames(
  short_sample_ids,
  metadata$sample_id
)

if (RUNTIME_CONFIG$debug_verbose) {
  message("\n=== Metadata and Color Scheme Initialization ===")
  # Metadata Loading
  message("\nMetadata Processing:")
  message(sprintf("  Source: %s", basename(metadata_path)))
  message(sprintf("  Rows: %d", nrow(metadata)))
  message(sprintf("  Columns: %d", ncol(metadata)))
  message("  Required columns present:")
  invisible(lapply(EXPERIMENT_CONFIG$COLUMN_ORDER, function(col) {
    message(sprintf("    %s: %s", col, col %in% colnames(metadata)))
  }))
  # Sample ID Mapping
  message("\nSample ID Processing:")
  message(sprintf("  Total IDs: %d", length(sample_ids)))
  message("  ID Mapping Examples (first 3):")
  head_ids <- head(names(sample_id_mapping), 3)
  invisible(lapply(head_ids, function(id) {
    message(sprintf("    %s -> %s", id, sample_id_mapping[id]))
  }))
  # Color Scheme
  message("\nColor Scheme Configuration:")
  message("  Fixed Colors:")
  message(sprintf("    Placeholder: %s", GENOME_TRACK_CONFIG$color_placeholder))
  message(sprintf("    Input: %s", GENOME_TRACK_CONFIG$color_input))
  message("\n  Category Colors:")
  message("    Antibody:")
  invisible(lapply(unique(metadata$antibody), function(ab) {
    message(sprintf("      %s: %s", ab, color_scheme$get_color("antibody", ab)))
  }))
  message("    Rescue Allele:")
  invisible(lapply(unique(metadata$rescue_allele), function(ra) {
    message(sprintf("      %s: %s", ra, color_scheme$get_color("rescue_allele", ra)))
  }))
  message("\n=== Initialization Complete ===")
}

# End message -------
message("Script completed succesfully...")
