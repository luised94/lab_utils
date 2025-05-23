#!/usr/bin/env Rscript
###############################################################################
# Plot bigwig files sequentially and in order
################################################################################
# PURPOSE: See all of the files in order
# USAGE: ./plot_genome_tracks_in_batches.R --experiment-id=<experiment-id> <options>
# DEPENDENCIES: GenomicRanges, rtracklayer
# OUTPUT: svg plots with comparisons for cpm/rpkm/raw and for shifted/raw/deduped for cpm counts
# AUTHOR: LEMR
# DATE: 2025-02-25
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
description <- "Plot genome tracks in batches"
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
    packages = required_packages,
    verbose = TRUE,
    skip_validation = args$skip_validation
)

################################################################################
# Load and Validate Experiment Configuration and Dependencies
################################################################################
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

# Setup logging if requested (independent)
#if (args$log_to_file) {
#    log_file <- setup_logging("chip_peak_analysis")
#    flog.appender(appender.file(log_file))
#    flog.threshold(INFO)
#}
################################################################################
# Setup directories, variables and metadata
################################################################################
# !! Add directories here.
required_directories <- c("fastq", "documentation", "coverage")
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

# Find bigwig files
bigwig_pattern <- sprintf("processed_.*_sequence_to_S288C_%s\\.bw$", EXPERIMENT_CONFIG$NORMALIZATION$active)
bigwig_files <- list.files(
    dirs$coverage,
    pattern = bigwig_pattern,
    full.names = TRUE
)

normalization_method <- sub(".*_([^_]+)\\.bw$", "\\1",
    basename(bigwig_files[1]))
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
        x = basename(bigwig_files)
    )
} else {
    # Extract sample IDs from fastq filenames
    sample_ids <- gsub(
        pattern = GENOME_TRACK_CONFIG$file_sample_id,
        replacement = "\\1",
        x = fastq_files
    )
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

# Create sample groups
sample_groups <- split(
    seq_len(nrow(metadata)),
    ceiling(seq_len(nrow(metadata)) / RUNTIME_CONFIG$process_samples_per_batch)
)

# Determine which groups to process
groups_to_process <- if (RUNTIME_CONFIG$process_single_file) {
    RUNTIME_CONFIG$process_batch
} else {
    seq_along(sample_groups)
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
    stop(sprintf("No reference genome files found matching pattern '%s' in: %s", 
                GENOME_TRACK_CONFIG$file_genome_pattern,
                GENOME_TRACK_CONFIG$file_genome_directory))
}

if (!file.exists(ref_genome_file)) {
    stop(sprintf("Reference genome file not accessible: %s", ref_genome_file[1]))
}
genome_data <- Biostrings::readDNAStringSet(ref_genome_file)

# Create chromosome range
chromosome_to_plot <- RUNTIME_CONFIG$process_chromosome
chromosome_width <- genome_data[chromosome_to_plot]@ranges@width
chromosome_roman <- paste0("chr", utils::as.roman(chromosome_to_plot))

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

# Calculate global track limits before entering for loop.
global_limits_result <- calculate_track_limits(
    bigwig_files = bigwig_files,
    genome_range = genome_range,
    padding_fraction = 0.1,
    verbose = RUNTIME_CONFIG$debug_verbose
)

if (!global_limits_result$success) {
    warning("Failed to calculate global y-limits: ", global_limits_result$error)
    global_limits <- GENOME_TRACK_CONFIG$track_ylim
} else {
    global_limits <- global_limits_result$data
}

# Process each group with all scaling modes
scaling_modes <- c("global", "local", "individual")

for (group_idx in groups_to_process) {
    if (RUNTIME_CONFIG$debug_verbose) {
        message("\n=== Processing Group ", group_idx, " ===")
        message(sprintf("Time: %s", format(Sys.time(), "%H:%M:%S")))
    }

    row_samples_to_visualize <- metadata[sample_groups[[group_idx]], ]

    if (RUNTIME_CONFIG$debug_verbose) {
        message("\nSample Selection:")
        message(sprintf("  Group Index: %d/%d", group_idx, length(groups_to_process)))
        message(sprintf("  Samples in group: %d", nrow(row_samples_to_visualize)))
        if (nrow(row_samples_to_visualize) > 0) {
            message("  Sample IDs:")
            invisible(lapply(row_samples_to_visualize$sample_id, 
                           function(id) message("    - ", id)))
        }
    }

    label_result <- create_track_labels(
        samples = row_samples_to_visualize,
        always_show = GENOME_TRACK_CONFIG$label_always_show,
        never_show = GENOME_TRACK_CONFIG$label_never_show,
        separator = GENOME_TRACK_CONFIG$label_separator,
        verbose = RUNTIME_CONFIG$debug_verbose
    )

    if (!label_result$success) {
        warning("Failed to create track labels for group ", group_idx, ": ", label_result$error)
        track_labels <- row_samples_to_visualize$short_name
        if (RUNTIME_CONFIG$debug_verbose) {
            message("  Using fallback labels (short_name)")
        }
    } else {
        track_labels <- label_result$data$labels
        if (RUNTIME_CONFIG$debug_verbose) {
            message("\nTrack Labels Created:")
            invisible(mapply(function(id, label) {
                message(sprintf("  %s -> %s", id, label))
            }, row_samples_to_visualize$sample_id, track_labels))
        }
    }
    # Initialize tracks list with chromosome axis
    tracks <- list(
        Gviz::GenomeAxisTrack(
            name = sprintf(GENOME_TRACK_CONFIG$format_genome_axis_track_name, chromosome_to_plot)
        )
    )

    for (i in seq_len(nrow(row_samples_to_visualize))) {
        sample_id <- row_samples_to_visualize$sample_id[i]
        current_antibody <- row_samples_to_visualize$antibody[i]
        bigwig_file_path <- bigwig_files[grepl(sample_id, bigwig_files)][1]
        if (RUNTIME_CONFIG$debug_verbose) {
            message("\nBigwig File Matching:")
            message(sprintf("  Sample ID: %s", sample_id))
            message("  Available files:")
            invisible(lapply(bigwig_files, function(f) message("    ", basename(f))))
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
            message(sprintf("  Track Label: %s", track_labels[i]))
            message(sprintf("  Bigwig file: %s", bigwig_file_path))
        }

        if (RUNTIME_CONFIG$debug_verbose) {
            message("\nTrack Creation Parameters:")
            message(sprintf("  Bigwig file: %s", 
                          if(is.na(bigwig_file_path)) "NOT FOUND" 
                          else basename(bigwig_file_path)))
            message("  Name arguments:")
            message(sprintf("    ID: %s", sample_id_mapping[sample_id]))
            message(sprintf("    Label: %s", track_labels[i]))
        }
        track_name_arguments <- c(
            sample_id_mapping[sample_id],  # Mapped ID
            track_labels[i]                # Generated label from create_track_labels
        )

        placeholder_name_arguments <- c(
            sample_id_mapping[sample_id],
            track_labels[i],
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

    plot_title <- sprintf(
        GENOME_TRACK_CONFIG$title_group_template,
        experiment_id,
        group_idx,
        chromosome_to_plot,
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

    for (mode in scaling_modes) {
        # Create filename for this mode
        plot_filename <- sprintf(
            GENOME_TRACK_CONFIG$filename_format_group_template,
            TIME_CONFIG$current_timestamp,
            experiment_id,
            group_idx,
            chromosome_to_plot,
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
            message(sprintf("    Experiment: %s", experiment_id))
            message(sprintf("    Chromosome: %s", chromosome_to_plot))
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
        if (mode == "global") {
            plot_config$ylim <- global_limits
        } else if (mode == "local") {
            # Calculate limits for this group's files
            group_files <- bigwig_files[grepl(
                paste(row_samples_to_visualize$sample_id, collapse = "|"),
                bigwig_files
            )]
            local_limits_result <- calculate_track_limits(
                bigwig_files = group_files,
                genome_range = genome_range,
                padding_fraction = .1,
                #padding_fraction = GENOME_TRACK_CONFIG$ylim_padding,
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
