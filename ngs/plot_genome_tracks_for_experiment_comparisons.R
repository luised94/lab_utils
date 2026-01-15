#!/usr/bin/env Rscript
################################################################################
# Plot bigwig files based EXPERIMENT_COMPARISONS configuration parameter
# Author: Luis | Date: 2025-05-02 | Version: 2.0.0
################################################################################
# PURPOSE: Plot specific comparisons based on quote expressions based on configuration_experiment_bmc of the experiment.
#
# USAGE:
#   ./core_scripts/plot_genome_tracks_for_experiment_comparisons.R --experiment-id=<experiment-id>
#
# DEPENDENCIES: 
#   ~/lab_utils/core_scripts/setup_bmc_experiment.R outputs
#   ~/lab_utils/core_scripts/{logging,script_control,file_operations}.R
#   required_packages
#
# OUTPUTS:
# - Svg or pdf files with genome tracks based on experiment comparisons written by the user
#
################################################################################
# Only run interactively
if(interactive()) {
  message("Running from repl... Loading functions.")
} else {
  stop("Run the script from the R repl in an interactive session.")
}

# Setup paths to configuration and root directory.
# Accounts for my reorganization of the repos.
ROOT_DIRECTORY <- system("git rev-parse --show-toplevel", intern = TRUE)
CORE_SCRIPTS_PATH <- file.path(ROOT_DIRECTORY, "core_scripts")
SCRIPT_CONFIGURATION_PATH <- file.path(
  CORE_SCRIPTS_PATH,
  "configuration_script_bmc.R"
)

#---------------------------------------
# Load user functions
#---------------------------------------
FUNCTION_FILENAME_TEMPLATE <- file.path(CORE_SCRIPTS_PATH, "functions_for_%s.R")
FUNCTION_FILENAMES <- c(
  "logging", "script_control",
  "file_operations", "bmc_config_validation",
  "metadata_processing", "genome_tracks"
)

for (function_filename in FUNCTION_FILENAMES) {
  function_filepath <- sprintf(
    FUNCTION_FILENAME_TEMPLATE,
    function_filename
  )
  normalized_path <- normalizePath(function_filepath)
  if (!file.exists(normalized_path)) {
    stop(sprintf("[FATAL] File with functions not found: %s", normalized_path))
  }
  source(normalized_path)
}

message("Loaded functions... Sourcing configuration.")

#---------------------------------------
# Source configuration for interactive session
#---------------------------------------
stopifnot(
  "Script configuration file does not exist. Please copy the template." =
    file.exists(SCRIPT_CONFIGURATION_PATH),
)
source(SCRIPT_CONFIGURATION_PATH)
message("Configuration file sourced... Checking configuration variables.")

# Ensure the variables expected in the script were //
# defined in the configuration file. //
# See template_interactive_script_configuration.R or //
# configuration_script_bmc.R //
required_configuration_variables <- c(
  "EXPERIMENT_ID", "EXPERIMENT_DIR",
  "CHROMOSOME_TO_PLOT", "OUTPUT_FORMAT",
  "OUTPUT_EXTENSION", "BIGWIG_PATTERN",
  "FASTQ_PATTERN", "SAMPLE_ID_CAPTURE_PATTERN",
  "ACCEPT_CONFIGURATION", "SKIP_PACKAGE_CHECKS"
)

is_missing_variable <- !sapply(required_configuration_variables, exists)
missing_variables <- required_configuration_variables[is_missing_variable]
if (length(missing_variables) > 0 ) {
  stop("Missing variable. Please define in 'script_configuration.R' file.",
       paste(missing_variables, collapse = ", "))
}

message("All variables defined in the configuration file...")


#-------------------------------------------------------------------------------
# Verify Required Libraries
#-------------------------------------------------------------------------------
# Add the packages that are used in the script.
required_packages <- c("rtracklayer", "GenomicRanges", "Gviz")
if (!is.character(required_packages) || length(required_packages) == 0) {
  stop("required_packages must be a non-empty character vector")
}

if (!SKIP_PACKAGE_CHECKS) {
  is_missing_package <- !sapply(
    X = required_packages,
    FUN = requireNamespace,
    quietly = TRUE
  )
  missing_packages <- required_packages[is_missing_package]
  if (length(missing_packages) > 0 ) {
    stop("Missing packages. Please install using renv:\n",
         paste(missing_packages, collapse = ", "))
  }
  SKIP_PACKAGE_CHECKS <- TRUE
}

message("All required packages available...")

#-------------------------------------------------------------------------------
# Setup experiment-specific configuration path, directories and file metadata
#-------------------------------------------------------------------------------
NUMBER_OF_EXPERIMENTS <- length(EXPERIMENT_DIR)
config_paths <- vector("character", length = NUMBER_OF_EXPERIMENTS)
metadata_paths <- vector("character", length = NUMBER_OF_EXPERIMENTS)
for (experiment_index in seq_len(NUMBER_OF_EXPERIMENTS)) {
  current_experiment_path <- EXPERIMENT_DIR[experiment_index]
  current_experiment_id <- EXPERIMENT_IDS[experiment_index]

  config_paths[experiment_index] <- file.path(
    current_experiment_path, "documentation",
    paste0(current_experiment_id, "_configuration_experiment_bmc.R")
  )

  metadata_paths[experiment_index] <- file.path(
    current_experiment_path, "documentation",
    paste0(current_experiment_id, "_sample_grid.csv")
  )
}

names(config_paths) <- EXPERIMENT_IDS
names(metadata_paths) <- EXPERIMENT_IDS

required_configuration_paths <- c(config_paths, metadata_paths)
is_missing_configuration_path <- !sapply(
  X = required_configuration_paths,
  FUN = file.exists
)
missing_configuration_paths <- required_configuration_paths[is_missing_configuration_path]

if ( length(missing_configuration_paths) > 0 ) {
  stop("Missing configuration paths. Please setup.\n",
       paste(missing_configuration_paths, collapse = ", "))
}

OUTPUT_DIR <- file.path(EXPERIMENT_DIR[1], "plots", "genome_tracks", "final_results")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
message("Configuration and metadata paths created. Loading metadata...")

#-------------------------------------------------------------------------------
# Setup directories, variables and metadata
#-------------------------------------------------------------------------------
# !! Add directories here.
required_directories <- c("fastq", "documentation", "coverage")
# Creates directories and stores them in dirs variable.
dirs <- setup_experiment_dirs(
    EXPERIMENT_DIR = EXPERIMENT_DIR,
    output_dir_name = "plots",
    required_input_dirs = required_directories
)

# Control output directory here.
dirs$output_dir <- file.path(dirs$plots, "genome_tracks", "experimental_comparisons")
dir.create(dirs$output_dir, recursive = TRUE, showWarnings = FALSE)

# Find fastq files and extract sample IDs
fastq_files <- list.files(
    path = dirs$fastq,
    pattern = GENOME_TRACK_CONFIG$file_pattern,
    full.names = FALSE
)

## Find bigwig files
#bigwig_pattern <- sprintf(
#    "processed_.*_sequence_to_S288C_blFiltered_%s\\.bw$",
#    EXPERIMENT_CONFIG$NORMALIZATION$active
#)

bigwig_files <- list.files(
    dirs$coverage,
    pattern = bigwig_pattern,
    full.names = TRUE
)

bigwig_basenames <- basename(bigwig_files)

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
        x = fastq_files
    )
    stopifnot(
        "Length of samples_ids is not lower than length of fastq files." =
        all(nchar(sample_ids) < nchar(fastq_files))
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

# Add after processing metadata but before track creation
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

#-------------------------------------------------------------------------------
# Setup genome and feature files
#-------------------------------------------------------------------------------
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
    stop(sprintf("No reference genome files found matching pattern '%s' in: %s",
                GENOME_TRACK_CONFIG$file_genome_pattern,
                GENOME_TRACK_CONFIG$file_genome_directory))
}

if (!file.exists(ref_genome_file)) {
    stop(sprintf("Reference genome file not accessible: %s", ref_genome_file[1]))
}
genome_data <- Biostrings::readDNAStringSet(ref_genome_file)

# Create chromosome range
chromosome_width <- genome_data[CHROMOSOME_TO_PLOT]@ranges@width
chromosome_roman <- paste0("chr", utils::as.roman(CHROMOSOME_TO_PLOT))

genome_range <- GenomicRanges::GRanges(
    seqnames = chromosome_roman,
    ranges = IRanges::IRanges(start = 1, end = chromosome_width),
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
                   GENOME_TRACK_CONFIG$file_feature_directory))
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
            dir.exists(GENOME_TRACK_CONFIG$file_genome_directory)),
        ".Feature Directory" = sprintf("%s (Exists: %s)",
            GENOME_TRACK_CONFIG$file_feature_directory,
            dir.exists(GENOME_TRACK_CONFIG$file_feature_directory)),

        "Genome Processing" = NULL,
        ".Search Pattern" = GENOME_TRACK_CONFIG$file_genome_pattern,
        ".Found File" = if(length(ref_genome_file) > 0) ref_genome_file else "None",
        ".File Accessible" = if(length(ref_genome_file) > 0) file.exists(ref_genome_file) else FALSE,

        "Chromosome Details" = NULL,
        ".Target" = CHROMOSOME_TO_PLOT,
        ".Roman Notation" = chromosome_roman,
        ".Width" = chromosome_width,

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

# Determine which comparisons to process
comparisons_to_process <- if (RUNTIME_CONFIG$process_single_comparison) {
    if (!RUNTIME_CONFIG$process_comparison %in% names(EXPERIMENT_CONFIG$COMPARISONS)) {
        cat(sprintf("Verify process_comparison in RUNTIME_CONFIG of configuration_experiment_bmc file for %s or verify override_configuration.R", EXPERIMENT_ID))
        stop("Debug comparison not found in EXPERIMENT_CONFIG$COMPARISONS")
    }
    RUNTIME_CONFIG$process_comparison
} else {
    names(EXPERIMENT_CONFIG$COMPARISONS)
}

if (RUNTIME_CONFIG$debug_verbose) {
    message("\nProcessing comparisons:")
    message(sprintf("- Total comparisons: %d", length(comparisons_to_process)))
    message("- Comparisons: ", paste(comparisons_to_process, collapse = ", "))
}

#-------------------------------------------------------------------------------
# Main script logic for plotting experiment comparisons of genomic tracks
#-------------------------------------------------------------------------------
scaling_modes <- c("local", "individual")
#scaling_modes <- "local"
for (comparison_name in comparisons_to_process) {
    # Metadata processing and sample selection
    comparison_expression <- EXPERIMENT_CONFIG$COMPARISONS[[comparison_name]]
    row_samples_to_visualize <- metadata[eval(comparison_expression,
        envir = metadata), ]

    # Try to find control for first available sample
    control_sample <- find_control_sample(
        row_samples_to_visualize[1, ],
        metadata,
        EXPERIMENT_CONFIG$CONTROL_FACTORS
    )

    # Combine the control sample to create the track labels.
    label_result <- create_track_labels(
        samples = rbind(row_samples_to_visualize, control_sample),
        always_show = GENOME_TRACK_CONFIG$label_always_show,
        never_show = GENOME_TRACK_CONFIG$label_never_show,
        separator = GENOME_TRACK_CONFIG$label_separator,
        verbose = RUNTIME_CONFIG$debug_verbose
    )

    if (!label_result$success) {
        warning("Failed to create track labels for comparison %s: %s", comparison_name, label_result$error)
        track_labels <- row_samples_to_visualize$short_name  # Fallback to sample_id
    } else {
        track_labels <- label_result$data$labels
    }

    if (nrow(row_samples_to_visualize) == 0) {
        warning(sprintf("No samples found for comparison: %s", comparison_name))
        next
    }

    # Initialize track list with genome axis
    tracks <- list(
        Gviz::GenomeAxisTrack(
           # name = paste("Chr ", CHROMOSOME_TO_PLOT, " Axis", sep = "")
            name = sprintf(GENOME_TRACK_CONFIG$format_genome_axis_track_name, CHROMOSOME_TO_PLOT)
        )
    )

    #------------------------------
    # Try to find control for first available sample
    control_sample <- find_control_sample(
        row_samples_to_visualize[1, ],
        metadata,
        EXPERIMENT_CONFIG$CONTROL_FACTORS
    )

    if(!is.null(control_sample)) {
        sample_id <- control_sample$sample_id
        current_antibody <- control_sample$antibody
        track_label <- track_labels[length(track_labels)]
        control_bigwig_file_path <- bigwig_files[grepl(control_sample$sample_id,
            bigwig_files)]

        track_name_arguments <- c(
            sample_id_mapping[sample_id],
            track_label
        )

        placeholder_name_arguments <- c(
            sample_id_mapping[sample_id],
            track_label,
            GENOME_TRACK_CONFIG$format_suffix
        )

        track_color <- if (current_antibody == "Input") {
            color_scheme$fixed$input
        } else {
            color_scheme$get_color("antibody", current_antibody)
        }

        control_track_creation_result <- create_control_track(
            bigwig_file_path = control_bigwig_file_path,
            track_format_name = GENOME_TRACK_CONFIG$format_sample_track_name,
            format_args = track_name_arguments,
            track_color = track_color,
            track_type = "h",
            genomic_range = genome_range,
            track_params = GENOME_TRACK_CONFIG$track_defaults_control,
            verbose = RUNTIME_CONFIG$debug_verbose
        )

        if (control_track_creation_result$success) {
            if (RUNTIME_CONFIG$debug_verbose) {
                message(sprintf("  Successfully created track for sample: %s", sample_id))
            }
            tracks[[length(tracks) + 1]] <- control_track_creation_result$data
        } else {
            # Report the error
            message(sprintf("\nTrack creation failed for sample %s:", sample_id))
            message(sprintf("  Error: %s", control_track_creation_result$error))
            message("  Creating placeholder track instead")
            if (RUNTIME_CONFIG$debug_verbose) {
                # Add context about the failure
                message("  Context:")
                message(sprintf("    Genomic Range: %s", 
                               if(exists("genome_range")) "Present" else "Missing"))
                message(sprintf("    Chromosome: %s", chromosome_roman))
                message(sprintf("    Width: %d", chromosome_width))
            }
            # Create placeholder with error handling
            placeholder_track_creation_result <- create_placeholder_track(
                sampling_rate = GENOME_TRACK_CONFIG$track_sampling_rate,
                chromosome_width = chromosome_width,
                track_color = GENOME_TRACK_CONFIG$color_placeholder,
                track_type = GENOME_TRACK_CONFIG$track_type,
                chromosome_name = chromosome_roman,
                placeholder_format_name = GENOME_TRACK_CONFIG$format_placeholder_track_name,
                format_args = placeholder_name_arguments,
                track_params = GENOME_TRACK_CONFIG$track_defaults_placeholder,
                verbose = RUNTIME_CONFIG$debug_verbose
            )
            if (!placeholder_track_creation_result$success) {
                stop(sprintf("Failed to create placeholder track: %s", placeholder_track_creation_result$error))
            } else {
                tracks[[length(tracks) + 1]] <- placeholder_track_creation_result$data
            }
        }
    }
    #------------------------------
    for (i in seq_len(nrow(row_samples_to_visualize))) {
        sample_id <- row_samples_to_visualize$sample_id[i]
        current_antibody <- row_samples_to_visualize$antibody[i]
        bigwig_file_path <- bigwig_files[grepl(sample_id, bigwig_files)]
        short_sample_id <- sample_id_mapping[sample_id]  # Mapped ID
        track_label <- track_labels[i]

        if (RUNTIME_CONFIG$debug_verbose) {
            message("\nBigwig File Matching:")
            message(sprintf("  Sample ID: %s", sample_id))
            message("  Available files:")
            invisible(lapply(bigwig_files[1:5], function(f) message("    ", basename(f))))
            message("  Pattern matches:")
            matches <- grepl(sample_id, bigwig_files)
            message(sprintf("    Matches found: %d", sum(matches)))
            if (sum(matches) > 0) {
                message("    Matching files:")
                invisible(lapply(bigwig_files[matches], function(f) 
                    message("      ", basename(f))))
            }
        }

        if (RUNTIME_CONFIG$debug_verbose) {
            message(sprintf("\nProcessing Sample %d/%d:", 
                          i, nrow(row_samples_to_visualize)))
            message(sprintf("  Sample ID: %s", sample_id))
            message(sprintf("  Antibody: %s", current_antibody))
            message(sprintf("  Track Label: %s", track_label))
            message(sprintf("  Bigwig file: %s", bigwig_file_path))
        }

        if (RUNTIME_CONFIG$debug_verbose) {
            message("\nTrack Creation Parameters:")
            message(sprintf("  Bigwig file: %s", 
                          if(is.na(bigwig_file_path)) "NOT FOUND" 
                          else basename(bigwig_file_path)))
            message("  Name arguments:")
            message(sprintf("    ID: %s", sample_id_mapping[sample_id]))
            message(sprintf("    Label: %s", track_label))
        }
        track_name_arguments <- c(
            short_sample_id,
            track_label
        )

        placeholder_name_arguments <- c(
            short_sample_id,
            track_label,
            GENOME_TRACK_CONFIG$format_suffix
        )

        track_color <- if (current_antibody == "Input") {
            color_scheme$fixed$input
        } else {
            color_scheme$get_color("antibody", current_antibody)
        }

        track_creation_result <- create_sample_track(
            bigwig_file_path = bigwig_file_path,
            track_format_name = GENOME_TRACK_CONFIG$format_sample_track_name,
            format_args = track_name_arguments,
            track_color = track_color,
            track_type = GENOME_TRACK_CONFIG$track_defaults_sample$type,
            genomic_range = genome_range,
            track_params = GENOME_TRACK_CONFIG$track_defaults_sample,
            verbose = RUNTIME_CONFIG$debug_verbose
        )

        if (track_creation_result$success) {
            if (RUNTIME_CONFIG$debug_verbose) {
                message(sprintf("  Successfully created track for sample: %s", sample_id))
            }
            tracks[[length(tracks) + 1]] <- track_creation_result$data
        } else {
            # Report the error
            message(sprintf("\nTrack creation failed for sample %s:", sample_id))
            message(sprintf("  Error: %s", track_creation_result$error))
            message("  Creating placeholder track instead")
            if (RUNTIME_CONFIG$debug_verbose) {
                # Add context about the failure
                message("  Context:")
                message(sprintf("    Genomic Range: %s", 
                               if(exists("genome_range")) "Present" else "Missing"))
                message(sprintf("    Chromosome: %s", chromosome_roman))
                message(sprintf("    Width: %d", chromosome_width))
            }
            # Create placeholder with error handling
            placeholder_track_creation_result <- create_placeholder_track(
                sampling_rate = GENOME_TRACK_CONFIG$track_sampling_rate,
                chromosome_width = chromosome_width,
                track_color = GENOME_TRACK_CONFIG$color_placeholder,
                type = GENOME_TRACK_CONFIG$track_type,
                chromosome_name = chromosome_roman,
                placeholder_format_name = GENOME_TRACK_CONFIG$format_placeholder_track_name,
                format_args = placeholder_name_arguments,
                track_params = GENOME_TRACK_CONFIG$track_defaults_placeholder,
                verbose = RUNTIME_CONFIG$debug_verbose
            )
            if (!placeholder_track_creation_result$success) {
                stop(sprintf("Failed to create placeholder track: %s", placeholder_track_creation_result$error))
            } else {
                tracks[[length(tracks) + 1]] <- placeholder_track_creation_result$data
            }
        }
    }
    #------------------------------

    # Add feature track if available
    if (exists("features")) {
        tracks[[length(tracks) + 1]] <- Gviz::AnnotationTrack(
            features,
            name = "Features",
            size = 0.5,
            background.title = "lightgray",
            fontcolor.title = "black",
            cex.title = 0.6
        )
    }

    base_comparison_name <- gsub("^comp_", "", comparison_name)  # Remove prefix
    sanitized_comparison_name <- gsub("_", ".", base_comparison_name)      # Convert underscores

    plot_title <- sprintf(
        GENOME_TRACK_CONFIG$title_comparison_template,
        EXPERIMENT_ID,
        sanitized_comparison_name,
        CHROMOSOME_TO_PLOT,
        nrow(row_samples_to_visualize),
        TIME_CONFIG$current_timestamp,
        normalization_method
    )

    plot_config <- create_track_plot_config(
        tracks = tracks,
        chromosome = chromosome_roman,
        to = chromosome_width,
        title = plot_title,
        visualization_params = GENOME_TRACK_CONFIG$plot_defaults,
        verbose = RUNTIME_CONFIG$debug_verbose
    )

    for (mode in scaling_modes[1]) {
        # Create filename for this mode
        plot_filename <- sprintf(
            GENOME_TRACK_CONFIG$filename_format_comparison_template,
            TIME_CONFIG$current_timestamp,
            EXPERIMENT_ID,
            sanitized_comparison_name,
            CHROMOSOME_TO_PLOT,
            mode
         )

         plot_file <- file.path(
             dirs$output_dir,
             plot_filename
         )

        if (RUNTIME_CONFIG$debug_verbose) {
            message(paste(rep("-", 80), collapse = ""))
            message(sprintf("\nScaling mode: %s", mode))
            message("\nPlot Generation Details:")
            message(sprintf("  Title Components:"))
            message(sprintf("    Experiment: %s", EXPERIMENT_ID))
            message(sprintf("    Chromosome: %s", CHROMOSOME_TO_PLOT))
            message(sprintf("    Sample Count: %d", nrow(row_samples_to_visualize)))
            message(sprintf("    Timestamp: %s", TIME_CONFIG$current_timestamp))
            message(sprintf("    Normalization: %s", normalization_method))
            message("\n  Output Configuration:")
            message(sprintf("    Plot Directory: %s", dirs$output_dir))
            message(sprintf("    Filename: %s", basename(plot_file)))
            message(sprintf("    Full Path: %s", plot_file))
            # Visual separator for readability in log
            message(paste(rep("-", 80), collapse = ""))
        }
        # Set y-limits based on mode
        if (mode == "local") {
            # Calculate limits for this group's files
            print(row_samples_to_visualize$sample_id)
            #group_files <- project_bigwig_files[grepl(
            #    paste(row_samples_to_visualize$sample_id, collapse = "|"),
            #    project_bigwig_files
            #)]
            pattern <- sprintf("processed_(%s)_sequence_to_S288C_blFiltered_CPM\\.bw$",
                               paste(row_samples_to_visualize$sample_id, collapse = "|"))
            group_files <- grep(pattern, bigwig_files, value = TRUE, perl = TRUE)
            if (RUNTIME_CONFIG$debug_verbose) {
                message("Selected sample IDs: ", paste(row_samples_to_visualize$sample_id, collapse=", "))
                message("Found files: ", paste(basename(group_files), collapse=", "))
                if (length(group_files) != nrow(row_samples_to_visualize)) {
                    warning(sprintf("Number of files (%d) doesn't match number of samples (%d)",
                            length(group_files), nrow(row_samples_to_visualize)))
                }
            }
            # INQ: Could cache this result as well and turn into independent script.
            local_limits_result <- calculate_track_limits(
                bigwig_files = group_files,
                genome_range = genome_range,
                padding_fraction = .1,
                verbose = RUNTIME_CONFIG$debug_verbose
            )
            plot_config$ylim <- if (local_limits_result$success) 
                local_limits_result$data else GENOME_TRACK_CONFIG$ylim_defaults
        } else {  # individual
            plot_config$ylim <- NULL  # Let tracks scale independently
        }

        if (RUNTIME_CONFIG$output_dry_run) {
            # Display only
            execute_track_plot(
                plot_config = plot_config,
                plot_params = GENOME_TRACK_CONFIG$plot_defaults,
                display_plot = TRUE,
                verbose = RUNTIME_CONFIG$debug_verbose
            )
        } else {
            # Save plot
            execute_track_plot(
                plot_config = plot_config,
                save_path = plot_file,
                save_params = list(
                    width = GENOME_TRACK_CONFIG$display_width,
                    height = GENOME_TRACK_CONFIG$display_height
                ),
                plot_params = GENOME_TRACK_CONFIG$plot_defaults,
                display_plot = FALSE,
                verbose = RUNTIME_CONFIG$debug_verbose
            )
        }
    }
}
