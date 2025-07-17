################################################################################
# BMC ChIP-seq Experiment Setup
# Author: Luis | Date: 2024-11-27 | Version: 2.0.0
################################################################################
#
# PURPOSE: Creates directory structure and metadata files for BMC ChIP-seq experiments
#
# USAGE:
# 1. Update experiment_id (format: YYMMDDBel, e.g., "241122Bel")
# 2. Source script to generate ~/data/[experiment_id]/ structure
#
# DEPENDENCIES: ~/lab_utils/core_scripts/bmc_config.R
#
# OUTPUTS:
# - Standard directory structure (peak/, fastq/, alignment/, bigwig/, plots/)
# - Sample metadata files (_sample_grid.csv, _bmc_table.tsv)
#
################################################################################
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
bigwig_pattern <- GENOME_TRACK_CONFIG$file_sample_id_from_bigwig
#bigwig_pattern <- sprintf(
#    "processed_.*_sequence_to_S288C_%s\\.bw$",
#    EXPERIMENT_CONFIG$NORMALIZATION$active
#)

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
    bigwig_basenames[1]
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
metadata_df <- load_and_process_experiment_metadata(
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
        antibody = unique(metadata_df$antibody),
        rescue_allele = unique(metadata_df$rescue_allele)
    )
)

# Shorten sample ids to minimal distinguishing set
short_sample_ids <- create_minimal_identifiers(
  metadata_df$sample_id,
  verbose = RUNTIME_CONFIG$debug_verbose
)

# Create mapping between full and short IDs
sample_id_mapping <- setNames(
  short_sample_ids,
  metadata_df$sample_id
)

if (RUNTIME_CONFIG$debug_verbose) {
  message("\n=== Metadata and Color Scheme Initialization ===")
  # Metadata Loading
  message("\nMetadata Processing:")
  message(sprintf("  Source: %s", basename(metadata_path)))
  message(sprintf("  Rows: %d", nrow(metadata_df)))
  message(sprintf("  Columns: %d", ncol(metadata_df)))
  message("  Required columns present:")
  invisible(lapply(EXPERIMENT_CONFIG$COLUMN_ORDER, function(col) {
    message(sprintf("    %s: %s", col, col %in% colnames(metadata_df)))
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
  invisible(lapply(unique(metadata_df$antibody), function(ab) {
    message(sprintf("      %s: %s", ab, color_scheme$get_color("antibody", ab)))
  }))
  message("    Rescue Allele:")
  invisible(lapply(unique(metadata_df$rescue_allele), function(ra) {
    message(sprintf("      %s: %s", ra, color_scheme$get_color("rescue_allele", ra)))
  }))
  message("\n=== Initialization Complete ===")
}
################################################################################
# Setup genome and feature files
################################################################################
stopifnot(
  "Genome directory not found" = dir.exists(GENOME_TRACK_CONFIG$file_genome_directory),
  "Feature directory not found" = dir.exists(GENOME_TRACK_CONFIG$file_feature_directory)
)

# Load reference genome
ref_genome_file <- list.files(
  path = GENOME_TRACK_CONFIG$file_genome_directory,
  pattern = GENOME_TRACK_CONFIG$file_genome_pattern,
  full.names = TRUE,
  recursive = TRUE
)[1]

if (length(ref_genome_file) == 0) {
  stop(sprintf(
        "No reference genome files found matching pattern '%s' in: %s",
        GENOME_TRACK_CONFIG$file_genome_pattern,
        GENOME_TRACK_CONFIG$file_genome_directory
  ))
}

if (!file.exists(ref_genome_file)) {
  stop(sprintf("Reference genome file not accessible: %s", ref_genome_file[1]))
}
genome_data <- Biostrings::readDNAStringSet(ref_genome_file)

# Create chromosome range
CHROMOSOME_TO_PLOT <- RUNTIME_CONFIG$process_chromosome
CHROMOSOME_WIDTH <- genome_data[CHROMOSOME_TO_PLOT]@ranges@width
CHROMOSOME_ROMAN <- paste0("chr", utils::as.roman(CHROMOSOME_TO_PLOT))

GENOME_RANGE_TO_LOAD <- GenomicRanges::GRanges(
    seqnames = CHROMOSOME_ROMAN,
    ranges = IRanges::IRanges(start = 1, end = CHROMOSOME_WIDTH),
    strand = "*"
)

# Load feature file (annotation)
feature_file <- list.files(
    path = GENOME_TRACK_CONFIG$file_feature_directory,
    pattern = GENOME_TRACK_CONFIG$file_feature_pattern,
    full.names = TRUE
)[1]

if (length(feature_file) == 0) {
  warning(sprintf("No feature files found matching pattern '%s' in: %s",
    GENOME_TRACK_CONFIG$file_feature_pattern,
    GENOME_TRACK_CONFIG$file_feature_directory
  ))
}

if (!is.null(feature_file)) {
  features <- rtracklayer::import(feature_file)
  # Convert to chrRoman format
  GenomeInfoDb::seqlevels(features) <- paste0(
    "chr",
    utils::as.roman(gsub("chr", "", GenomeInfoDb::seqlevels(features)))
  )
}

if (RUNTIME_CONFIG$debug_verbose) {
  debug_info <- list(
    "title" = "Genome and Feature Loading Status",
    "Directory Validation" = NULL,
    ".Genome Directory" = sprintf("%s (Exists: %s)",
      GENOME_TRACK_CONFIG$file_genome_directory,
      dir.exists(GENOME_TRACK_CONFIG$file_genome_directory
      )),
    ".Feature Directory" = sprintf("%s (Exists: %s)",
      GENOME_TRACK_CONFIG$file_feature_directory,
      dir.exists(GENOME_TRACK_CONFIG$file_feature_directory
      )),
    "Genome Processing" = NULL,
    ".Search Pattern" = GENOME_TRACK_CONFIG$file_genome_pattern,
    ".Found File" = if(length(ref_genome_file) > 0) ref_genome_file else "None",
    ".File Accessible" = if(length(ref_genome_file) > 0) file.exists(ref_genome_file) else FALSE,
    "Chromosome Details" = NULL,
    ".Target" = CHROMOSOME_TO_PLOT,
    ".Roman Notation" = CHROMOSOME_ROMAN,
    ".Width" = CHROMOSOME_WIDTH,
    "Feature Processing" = NULL,
    ".Feature Search Pattern" = GENOME_TRACK_CONFIG$file_feature_pattern,
    ".Feature Found File" = if(length(feature_file) > 0) feature_file else "None",
    ".Feature File Accessible" = if(!is.null(feature_file) && length(feature_file) > 0)
            file.exists(feature_file) else FALSE
    )
    # Add feature details if loaded
    if (exists("features") && !is.null(features)) {
        debug_info[[".Feature Count"]] <- length(features)
        debug_info[[".Sequence Levels"]] <- paste(GenomeInfoDb::seqlevels(features), collapse = ", ")
    }
    print_debug_info(debug_info)
}


################################################################################
# MAIN
################################################################################
# Filter all of the bigwig files to CPM only.
is_CPM_bw_file <- grepl("CPM\\.bw$", bigwig_files)
bigwig_files_subset <- bigwig_files[is_CPM_bw_file]

#bigwig_file_count_per_sample <- unlist(lapply(seq_len(nrow(metadata_df)) function(row_idx){
#  is_bw_of_sample_id <- grepl(metadata_df[row_idx, "sample_id"], bigwig_files_subset)
#  return(sum(is_bw_of_sample_id))
#}))

# Initialize track list with genome axis
track_container <- list(
  Gviz::GenomeAxisTrack(
  name = sprintf(
    GENOME_TRACK_CONFIG$format_genome_axis_track_name,
    CHROMOSOME_TO_PLOT
    )
  )
)
#if (exists("features")) {
#  track_container[[length(track_container) + 1]] <- Gviz::AnnotationTrack(
#  features,
#  name = "Features",
#  size = 0.5,
#  background.title = "lightgray",
#  fontcolor.title = "black",
#  cex.title = 0.6
#  )
#}

plot_prefix <- "blacklist_processing_effect"
MAX_ROW <- nrow(metadata_df)
#ROWS_IDX_TO_PLOT <- seq_len(MAX_ROW)[1:3]
ROWS_IDX_TO_PLOT <- seq_len(MAX_ROW)
for (row_idx in ROWS_IDX_TO_PLOT) {
  message("--- For loop for metadata ---")
  message(
    sprintf("  Processing row: %s / %s ",
    row_idx, MAX_ROW)
  )
  current_row_df <- metadata_df[row_idx, ]
  current_sample_id <- current_row_df$sample_id
  is_bw_of_sample_id <- grepl(current_sample_id, bigwig_files_subset)
  current_bigwig_files_subset <- bigwig_files_subset[is_bw_of_sample_id]
  track_color <- color_scheme$get_color("antibody", current_row_df[, "antibody"])

  message(
    sprintf("  Current sample id: %s\n  Amount of bigwig files to plot: %s\n Track color: %s",
      current_sample_id, length(current_bigwig_files_subset), track_color)
  )
  #message("  Current row:")
  #print(current_row_df[, c(EXPERIMENT_CONFIG$COLUMN_ORDER, "sample_id")])

  # Grab the first track to reset the list container.
  track_addition_count <- 0
  track_container <- list(track_container[[1]])
  for (bigwig_file_path in current_bigwig_files_subset) {
    message("  --- For loop for bigwig file ---")
    parts_of_bigwig_file_path <- strsplit(x = bigwig_file_path, split = "_", fixed = TRUE)
    bigwig_type <- parts_of_bigwig_file_path[[1]][length(parts_of_bigwig_file_path[[1]])-1]
    track_name_arguments <- c(
      sample_id_mapping[current_sample_id],
      as.character(current_row_df[, "antibody"]),
      bigwig_type
    )
    track_name <- paste(track_name_arguments, collapse = ".")

    message("  Current bigwig file: ", bigwig_file_path)
    message("  Sample id mapping: ", sample_id_mapping[current_sample_id])
    message("  Row idx: ", as.numeric(row_idx))
    #message("  Parts of bigwig_file_path: ")
    #print(parts_of_bigwig_file_path)
    message("  Length of parts vector: ", as.character(length(parts_of_bigwig_file_path[[1]])))
    message("  Next to last part: ", bigwig_type)
    message("  Track name: ", track_name)
    message("  ~~~~~~~~~~~~~~~~~~~~~~~~~")

    current_bigwig_gr <- rtracklayer::import(bigwig_file_path, which = GENOME_RANGE_TO_LOAD)
    track_data <- Gviz::DataTrack(
      range = current_bigwig_gr,
      name = track_name,
      col = track_color,
      type = "h",
      size = 1.2,
      showaxis = TRUE,
      showtitle = TRUE,
      background.title = "white",
      fontcolor.title = "black",
      col.border.title = "#e0e0e0",
      cex.title = 0.6,
      fontface = 1,
      title.width = 1.2
    )
    #track_container <- append(
    #  x = track_container,
    #  values = track_data,
    #  after = 1 + track_addition_count + 1
    #)
    track_container[[length(track_container) + 1]] <- track_data

  }

  current_row_values_for_name <- lapply(
    EXPERIMENT_CONFIG$COLUMN_ORDER,
    function(column_name){
      as.character(current_row_df[, column_name])
    }
  )
  plot_name <- paste0(
    paste(plot_prefix, current_sample_id, row_idx, sep = "_"),
    "_",
    paste(current_row_values_for_name, collapse = "."),
    ".svg"
  )
  plot_title <- paste(plot_prefix, current_sample_id, row_idx, sep = "_")
  plot_output_path <- file.path(dirs$output_dir, plot_name)
  message("    Plot name: ", plot_name)
  message("    Plot title: ", plot_title)
  message("    Plot output path: ", plot_output_path)
  # Add feature track if available
  if (exists("features")) {
    track_container[[length(track_container) + 1]] <- Gviz::AnnotationTrack(
      features,
      name = "Features",
      size = 0.5,
      background.title = "lightgray",
      fontcolor.title = "black",
      cex.title = 0.6
    )
  }

  if (file.exists(plot_output_path)) {
    message("Plot output already exists. Skipping...")
  } else {
    svglite::svglite(
      filename = plot_output_path,
      width = 10,
      height = 8,
      bg = "white"
    )
    Gviz::plotTracks(
      trackList = track_container,
      chromosome = CHROMOSOME_ROMAN,
      from = GENOME_RANGE_TO_LOAD@ranges@start,
      to = GENOME_RANGE_TO_LOAD@ranges@width,
      margin = 15,
      innerMargin = 5,
      spacing = 10,
      main = plot_title,
      col.axis = "black",
      cex.axis = 0.8,
      cex.main = 0.7,
      fontface.main = 1,
      background.panel = "transparent"
    )
    dev.off()
  message("   Plot saved...")
  }
  message("\n")
  message("=========================")
}
# End message -------
message("Script completed succesfully...")
