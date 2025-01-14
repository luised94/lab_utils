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
    ),
    list(
        path = "~/lab_utils/core_scripts/functions_for_genome_tracks.R",
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

if (length(fastq_files) == 0) {
    stop("No fastq files found in specified directory")
}

# Extract sample IDs from fastq filenames
sample_ids <- gsub(
    pattern = GENOME_TRACK_CONFIG$file_sample_id,
    replacement = "\\1",
    x = fastq_files
)

# Find bigwig files
bigwig_pattern <- sprintf("_%s\\.bw$", EXPERIMENT_CONFIG$NORMALIZATION$active)
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
    message("\nColor scheme initialized:")
    message("Fixed colors:")
    print(color_scheme$fixed)
    message("Category colors:")
    print(color_scheme$categories)
}

# Create sample groups
sample_groups <- split(
    seq_len(nrow(metadata)),
    ceiling(seq_len(nrow(metadata)) / RUNTIME_CONFIG$process_samples_per_batch)
)

# Determine which groups to process
groups_to_process <- if (RUNTIME_CONFIG$process_single_file) {
    RUNTIME_CONFIG$process_group
} else {
    seq_along(sample_groups)
}

################################################################################
# Setup genome and feature files
################################################################################
# Load reference genome
ref_genome_file <- list.files(
    GENOME_TRACK_CONFIG$file_genome_directory,
    pattern = GENOME_TRACK_CONFIG$file_genome_pattern,
    full.names = TRUE,
    recursive = TRUE
)[1]
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
    GENOME_TRACK_CONFIG$file_feature_directory,
    pattern = GENOME_TRACK_CONFIG$file_feature_pattern,
    full.names = TRUE
)[1]

if (!is.null(feature_file)) {
    features <- rtracklayer::import(feature_file)
    # Convert to chrRoman format
    GenomeInfoDb::seqlevels(features) <- paste0(
        "chr",
        utils::as.roman(gsub("chr", "", GenomeInfoDb::seqlevels(features)))
    )
}

limits_result <- calculate_track_limits(
    bigwig_files = bigwig_files,
    genome_range = genome_range,
    padding_fraction = 0.1,
    verbose = RUNTIME_CONFIG$debug_verbose
)

if (!limits_result$success) {
    warning("Failed to calculate y-limits: ", limits_result$error)
    y_limits <- c(0, 1000)  # Default fallback
} else {
    y_limits <- limits_result$data
}

#for (group_idx in groups_to_process) {
#
#    row_samples_to_visualize <- metadata[sample_groups[[group_idx]], ]
#
#    label_result <- create_track_labels(
#        samples = row_samples_to_visualize,
#        always_show = "antibody",
#        never_show = c("sample_id", "full_name", "short_name", "X__cf_genotype"),
#        separator = "-",
#        verbose = TRUE
#    )
#
#    if (!label_result$success) {
#        warning("Failed to create track labels for comparison %s: %s", comparison_name, label_result$error)
#        track_labels <- row_samples_to_visualize$short_name  # Fallback to sample_id
#    } else {
#        track_labels <- label_result$data$labels
#    }
#
#    if (nrow(row_samples_to_visualize) == 0) {
#        warning(sprintf("No samples found for comparison: %s", comparison_name))
#        next
#    }
#
#    if (RUNTIME_CONFIG$debug_verbose) {
#        message(sprintf("Found %d samples for comparison", nrow(row_samples_to_visualize)))
#    }
#
#    # Initialize tracks list with chromosome axis
#    tracks <- list(
#        Gviz::GenomeAxisTrack(
#            name = sprintf(track_axis_name_template, chromosome_to_plot)
#        )
#    )
#
#    for (i in seq_len(nrow(row_samples_to_visualize))) {
#        sample_id <- row_samples_to_visualize$sample_id[i]
#        current_antibody <- row_samples_to_visualize$antibody[i]
#        # Find matching bigwig file
#        bigwig_file_path <- bigwig_files[grepl(sample_id, bigwig_files)][1]
#
#        track_name <- sprintf(
#            GENOME_TRACK_CONFIG$format_track,
#            sample_id_mapping[sample_id],
#            row_samples_to_visualize$short_name[i],
#            row_samples_to_visualize$antibody[i]
#        )
#        track_name_arguments <- c(
#            sample_id_mapping[sample_id],
#            row_samples_to_visualize$short_name[i],
#            row_samples_to_visualize$antibody[i]
#        )
#
#        placeholder_name <- sprintf(
#            "%s %s",
#            track_name,
#            GENOME_TRACK_CONFIG$format_placeholder
#        )
#        track_color <- if (current_antibody == "Input") {
#            color_scheme$fixed$input
#        } else {
#            color_scheme$get_color("antibody", current_antibody)
#        }
#        track_creation_result <- create_sample_track(
#            bigwig_file_path = bigwig_file_path,
#            track_format_name = ,
#            format_args = track_name_arguments,
#            track_color = track_color,
#            track_type = GENOME_TRACK_CONFIG$track_type,
#            genomic_range = genomic_range
#
#        )
#        if (track_creation_result$success) {
#            tracks[[length(tracks) + 1]] <- track_creation_result$data
#        } else {
#            track_name_arguments <- c(
#                sample_id_mapping[sample_id],
#                GENOME_TRACK_CONFIG$format_suffix
#            )
#            placeholder_track_creation_result <- create_placeholder_track(
#                sampling_rate = 100,
#                chromosome_width = chromosome_width,
#                track_color = track_color,
#                type = GENOME_TRACK_CONFIG$track_type,
#                chromosome_name = chromosome_roman,
#                placeholder_format_name = GENOME_TRACK_CONFIG$format_placeholder_track_name,
#                format_args = track_name_arguments
#            )
#            if(placeholder_track_creation_result$success) {
                 #if (RUNTIME_CONFIG$debug_verbose) {
                 #    message("Created placeholder for sample: ", sample_id)
                 #}
#            tracks[[length(tracks) + 1]] <- placeholder_track_creation_result$data
#           } else {
#                stop(sprintf("Placeholder track creation failed for sample id %s.", sample_id))
#           }
#        }
#    }
#    # Add feature track if available
#    if (exists("features")) {
#        tracks[[length(tracks) + 1]] <- Gviz::AnnotationTrack(
#            features,
#            name = "Features"
#        )
#    }
#    plot_title <- sprintf(
#        GENOME_TRACK_CONFIG$title_dev_template,
#        experiment_id,
#        chromosome_to_plot,
#        nrow(row_samples_to_visualize),
#        TIMESTAMPS$full,
#        normalization_method
#    )
#    #plot_config <- create_track_plot_config(
#    #    tracks = tracks,
#    #    chromosome = chromosome_roman,
#    #    to = chromosome_width,
#    #    ylim = y_limits,
#    #    title = plot_title,
#    #    visualization_params = viz_params
#    #)
#    #execute_track_plot(
#    #    plot_config = plot_config,
#    #    save_path = NULL,
#    #    save_params = list()
#
#    #)
#    #if (RUNTIME_CONFIG$output_save_plots) {
#   # 
#   #     plot_filename <- sprintf(
#   #             "%s_%s_chr%s_n%d_group%d.svg",
#   #             TIMESTAMPS$full,
#   #             experiment_id,
#   #             chromosome_to_plot,
#   #             nrow(row_samples_to_visualize),
#   #             group_idx
#   #     )
#   #     plot_file <- file.path(
#   #         plots_dir,
#   #         plot_filename
#   #     )
#        #if (RUNTIME_CONFIG$debug_verbose) {
#        #    message(sprintf(
#        #        "Saving plot to: %s",
#        #        basename(plot_file)
#        #    ))
#        #}
#    #    execute_track_plot(
#    #        plot_config,
#    #        save_path = plot_file,
#    #        save_params = list(
#    #            width = GENOME_TRACK_CONFIG$display_width,
#    #            height = GENOME_TRACK_CONFIG$display_height
#    #        )
#    #    )
#    #}
    ## Interactive viewing options
    #if (RUNTIME_CONFIG$debug_interactive) {
    #    user_input <- readline(
    #        prompt = GENOME_TRACK_CONFIG$interactive_prompt
    #    )
    #    if (user_input == "q") break
    #    if (user_input == "s") RUNTIME_CONFIG$output_save_plots <- FALSE
    #} else {
    #    Sys.sleep(RUNTIME_CONFIG$output_display_time)  # Pause between plots
    #}
#}
