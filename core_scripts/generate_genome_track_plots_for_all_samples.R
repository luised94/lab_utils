#!/usr/bin/env Rscript
#REPLACE
source(file.path(Sys.getenv("HOME"), "lab_utils", "core_scripts", "functions_for_script_control.R"))

# Parse arguments and validate configurations
args <- parse_args(commandArgs(trailingOnly = TRUE))
experiment_id <- args[["experiment-id"]]
source(file.path("~/data", experiment_id, "documentation", 
                paste0(experiment_id, "_bmc_config.R")))
validate_configs(c("RUNTIME_CONFIG", "EXPERIMENT_CONFIG"))
################################################################################
# Load Required Libraries
################################################################################
required_packages <- c("rtracklayer", "GenomicRanges", "Gviz")
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' is missing", pkg))
    }
}

################################################################################
# Load and Validate Experiment Configuration and Dependencies
################################################################################
#source("~/lab_utils/core_scripts/all_functions.R")
source("~/lab_utils/core_scripts/functions_for_genome_tracks.R")
source("~/lab_utils/core_scripts/functions_for_metadata_processing.R")
# Load metadata and files
#-----------------------------------------------------------------------------
#experiment_id <- "241007Bel"
experiment_id <- "100303Bel"
base_dir <- file.path(Sys.getenv("HOME"), "data", experiment_id)
plots_dir <- file.path(base_dir, "plots", "genome_tracks", "overview")
metadata_path <- file.path(base_dir, "documentation",
                          paste0(experiment_id, "_sample_grid.csv"))
config_path <- file.path(base_dir, "documentation", paste0(experiment_id, "_bmc_config.R"))
source(config_path)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Find fastq files and extract sample IDs
fastq_files <- list.files(
    path = file.path(base_dir, "fastq"),
    pattern = "consolidated_.*_sequence\\.fastq$",
    full.names = FALSE
)

if (length(fastq_files) == 0) {
    stop("No fastq files found in specified directory")
}

# Extract sample IDs from fastq filenames
sample_ids <- gsub(
    pattern = "consolidated_([0-9]{5,6})_sequence\\.fastq",
    replacement = "\\1",
    x = fastq_files
)

# 3. Load and process metadata
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)

# 4. Enforce factor levels from config
for (col_name in names(EXPERIMENT_CONFIG$CATEGORIES)) {
    if (col_name %in% colnames(metadata)) {
        metadata[[col_name]] <- factor(
            metadata[[col_name]],
            levels = EXPERIMENT_CONFIG$CATEGORIES[[col_name]],
            ordered = TRUE
        )
    }
}

# 5. Sort metadata using config column order
sorted_metadata <- metadata[do.call(
    order,
    metadata[EXPERIMENT_CONFIG$COLUMN_ORDER]
), ]

# 6. Add sample IDs to metadata
sorted_metadata$sample_id <- sample_ids

if (RUNTIME_CONFIG$debug_verbose) {
    message("Metadata processing summary:")
    message(sprintf("Found %d fastq files", length(fastq_files)))
    message(sprintf("Processed %d metadata rows", nrow(sorted_metadata)))
    message("Columns with enforced factors:")
    print(names(EXPERIMENT_CONFIG$CATEGORIES))
}

# Create color scheme
color_scheme <- create_color_scheme(
    config = list(
        placeholder = GENOME_TRACK_CONFIG$color_placeholder,
        input = GENOME_TRACK_CONFIG$color_input
    ),
    categories = list(
        antibody = unique(sorted_metadata$antibody),
        rescue_allele = unique(sorted_metadata$rescue_allele)
    )
)

if (RUNTIME_CONFIG$debug_verbose) {
    message("\nColor scheme initialized:")
    message("Fixed colors:")
    print(color_scheme$fixed)
    message("Category colors:")
    print(color_scheme$categories)
}


# Load reference genome
ref_genome_file <- list.files(
    file.path(Sys.getenv("HOME"), "data", "REFGENS"),
    pattern = "S288C_refgenome.fna",
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

# Find bigwig files
bigwig_pattern <- sprintf("_%s\\.bw$", EXPERIMENT_CONFIG$NORMALIZATION$active)
bigwig_files <- list.files(
    file.path(base_dir, "coverage"),
    pattern = bigwig_pattern,
    full.names = TRUE
)
normalization_method <- sub(".*_([^_]+)\\.bw$", "\\1",
                          basename(bigwig_files[1]))

# Load feature file (annotation)
feature_file <- list.files(
    file.path(Sys.getenv("HOME"), "data", "feature_files"),
    pattern = "eaton_peaks",
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

# Process samples in groups
#-----------------------------------------------------------------------------
# Create sample groups
sample_groups <- split(
    seq_len(nrow(sorted_metadata)),
    ceiling(seq_len(nrow(sorted_metadata)) / RUNTIME_CONFIG$process_samples_per_group)
)

# Determine which groups to process
groups_to_process <- if (RUNTIME_CONFIG$process_single_file) {
    RUNTIME_CONFIG$process_group
} else {
    seq_along(sample_groups)
}

# Calculate global y-limits for all plots (before the plotting loop)
#-----------------------------------------------------------------------------
if (RUNTIME_CONFIG$debug_verbose) {
    message("\nCalculating global range for all tracks...")
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
# Add after processing metadata but before track creation
short_sample_ids <- create_minimal_identifiers(sorted_metadata$sample_id, verbose = RUNTIME_CONFIG$debug_verbose)
# Create mapping between full and short IDs
sample_id_mapping <- setNames(short_sample_ids, sorted_metadata$sample_id)

# Process each group
for (group_idx in groups_to_process) {
    if (RUNTIME_CONFIG$debug_verbose) {
        message("\nProcessing group ", group_idx)
    }

    #===============================================================================
    # 1. Metadata Processing and Sample Selection
    #===============================================================================

    # Get current group's samples
    row_samples_to_visualize <- sorted_metadata[sample_groups[[group_idx]], ]

    label_result <- create_track_labels(
        samples = row_samples_to_visualize,
        always_show = "antibody",
        never_show = c("sample_id", "full_name", "short_name", "X__cf_genotype"),
        separator = "-",
        verbose = TRUE
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

    if (RUNTIME_CONFIG$debug_verbose) {
        message(sprintf("Found %d samples for comparison", nrow(row_samples_to_visualize)))
    }

    # Initialize tracks list with chromosome axis
    tracks <- list(
        Gviz::GenomeAxisTrack(
            name = sprintf(track_axis_name_template, chromosome_to_plot)
        )
    )

    # Add tracks for each sample
    for (i in seq_len(nrow(row_samples_to_visualize))) {
        sample_id <- row_samples_to_visualize$sample_id[i]
        current_antibody <- row_samples_to_visualize$antibody[i]

        # Find matching bigwig file
        bigwig_file_path <- bigwig_files[grepl(sample_id, bigwig_files)][1]

        track_name <- sprintf(
            GENOME_TRACK_CONFIG$format_track,
            sample_id_mapping[sample_id],
            row_samples_to_visualize$short_name[i],
            row_samples_to_visualize$antibody[i]
        )

        placeholder_name <- sprintf(
            "%s %s",
            track_name,
            GENOME_TRACK_CONFIG$format_placeholder
        )
        track_color <- if (current_antibody == "Input") {
            color_scheme$fixed$input
        } else {
            color_scheme$get_color("antibody", current_antibody)
        }
        #track_creation_result <- create_sample_track(
        #    bigwig_file_path,
        #    track_format_name = GENOME_TRACK_CONFIG$format_control,
        #    format_args = c(sample_id_mapping[sample_id], short_name, current_antibody)
        #    track_color = color_scheme,
        #    track_type = track_type,
        #    genomic_range = genomic_range

        #)
        #if(track_creation_result$success) {
        #    tracks[[length(tracks) + 1]] <- track_result$data
        #} else {
        #    # Create placeholder
        #    tracks[[length(tracks) + 1]] <- create_placeholder_track(
        #        sampling_rate = 100,
        #        chromosome_width = NULL,
        #        track_color,
        #        type = "l",
        #        chromosome_name,
        #        placeholder_format_name,
        #        format_args
        #   )
        #}
        if (!is.na(bigwig_file_path) && file.exists(bigwig_file_path)) {
            if (RUNTIME_CONFIG$debug_verbose) {
                message("Adding track for sample: ", sample_id)
            }

            # Update track color based on antibody
            track_color <- if (current_antibody == "Input") {
                color_scheme$fixed$input
            } else {
                color_scheme$get_color("antibody", current_antibody)
            }
            track_data <- rtracklayer::import(bigwig_file_path)

            tracks[[length(tracks) + 1]] <- Gviz::DataTrack(
                track_data,
                name = track_name,
                type = "l",
                col = track_color
            )
        } else {
            if (RUNTIME_CONFIG$debug_verbose) {
                message("Creating placeholder for sample: ", sample_id)
            }

            # Calculate reasonable number of points (e.g., one point per 100bp)
            sampling_rate <- 100  # Adjust this value based on your needs
            num_points <- ceiling(chromosome_width / sampling_rate)

            # Create proper placeholder track with genome coordinates
            empty_ranges <- GenomicRanges::GRanges(
                seqnames = chromosome_roman,
                ranges = IRanges::IRanges(
                    # Create evenly spaced points across chromosome
                    start = seq(1, chromosome_width, length.out = num_points),
                    width = 1
                ),
                score = rep(0, num_points),  # One score per position
                strand = "*"
            )

            empty_track <- Gviz::DataTrack(
                empty_ranges,
                name = placeholder_name,
                type = "l",
                col = color_scheme$fixed$placeholder,
                chromosome = chromosome_roman
            )
            tracks[[length(tracks) + 1]] <- empty_track
        }
    }

    # Add feature track if available
    if (exists("features")) {
        tracks[[length(tracks) + 1]] <- Gviz::AnnotationTrack(
            features,
            name = "Features"
        )
    }

   # viz_params <- if (GENOME_TRACK_CONFIG$use_custom_visualization) {
   #     list(
   #         fontcolor = GENOME_TRACK_CONFIG$track_fontcolor,
   #         background.title = GENOME_TRACK_CONFIG$track_background
   #         # Other parameters...
   #     )
   # } else {
   #     list()  # Empty list will use all defaults
   # }
    plot_title <- sprintf(
        GENOME_TRACK_CONFIG$title_dev_template,
        experiment_id,
        chromosome_to_plot,
        nrow(row_samples_to_visualize),
        TIMESTAMPS$full,
        normalization_method
    )
    #plot_config <- create_track_plot_config(
    #    tracks = tracks,
    #    chromosome = chromosome_roman,
    #    to = chromosome_width,
    #    ylim = y_limits,
    #    title = plot_title,
    #    visualization_params = viz_params
    #)
    #execute_track_plot(
    #    plot_config = plot_config,
    #    save_path = NULL,
    #    save_params = list()

    #)
    #if (RUNTIME_CONFIG$output_save_plots) {
   # 
   #     plot_filename <- sprintf(
   #             "%s_%s_chr%s_n%d_group%d.svg",
   #             TIMESTAMPS$full,
   #             experiment_id,
   #             chromosome_to_plot,
   #             nrow(row_samples_to_visualize),
   #             group_idx
   #     )
   #     plot_file <- file.path(
   #         plots_dir,
   #         plot_filename
   #     )
    #    execute_track_plot(
    #        plot_config,
    #        save_path = plot_file,
    #        save_params = list(
    #            width = GENOME_TRACK_CONFIG$display_width,
    #            height = GENOME_TRACK_CONFIG$display_height
    #        )
    #    )
    #}
    Gviz::plotTracks(
        trackList = tracks,
        chromosome = chromosome_roman,
        from = 1,
        to = chromosome_width,
        ylim = y_limits
        #main = plot_title,
        ## Track name appearance
        #fontcolor = "black",           # Track name text color
        #background.title = "white",    # Track name background
        #col.border.title = "#E0E0E0",  # Light gray border around track names

        ## Other visualization parameters
        #cex.main = 1,
        #margin = 15,
        #innerMargin = 5,

        ## Axis appearance
        #col.axis = "black",            # Axis text color
        #cex.axis = 0.8,               # Axis text size

        ## Title panel
        #col.title = "black",          # Title text color
        #fontface.title = 2            # Bold title
    )

    # Save plot if needed
    if (RUNTIME_CONFIG$output_save_plots) {
        plot_filename <- sprintf(
                "%s_%s_chr%s_n%d_group%d.svg",
                TIMESTAMPS$full,
                experiment_id,
                chromosome_to_plot,
                nrow(row_samples_to_visualize),
                group_idx
        )
        plot_file <- file.path(
            plots_dir,
            plot_filename
        )

        if (RUNTIME_CONFIG$debug_verbose) {
            message(sprintf(
                "Saving plot to: %s",
                basename(plot_file)
            ))
        }

        svg(plot_file, width = GENOME_TRACK_CONFIG$width, height = GENOME_TRACK_CONFIG$display_height)
        Gviz::plotTracks(
            trackList = tracks,
            chromosome = chromosome_roman,
            from = 1,
            to = chromosome_width,
            ylim = y_limits,
            main = plot_title,
            # Track name appearance
            fontcolor = "black",           # Track name text color
            background.title = "white",    # Track name background
            col.border.title = "#E0E0E0",  # Light gray border around track names

            # Other visualization parameters
            cex.main = 1,
            margin = 15,
            innerMargin = 5,

            # Axis appearance
            col.axis = "black",            # Axis text color
            cex.axis = 0.8,               # Axis text size

            # Title panel
            col.title = "black",          # Title text color
            fontface.title = 2            # Bold title
        )
        dev.off()
    }

    # Interactive viewing options
    if (RUNTIME_CONFIG$debug_interactive) {
        user_input <- readline(
            prompt = GENOME_TRACK_CONFIG$interactive_prompt
        )
        if (user_input == "q") break
        if (user_input == "s") RUNTIME_CONFIG$output_save_plots <- FALSE
    } else {
        Sys.sleep(RUNTIME_CONFIG$output_display_time)  # Pause between plots
    }
}

message("Processing complete")
