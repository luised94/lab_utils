#!/usr/bin/env Rscript
# Configuration
#-----------------------------------------------------------------------------
DEBUG_CONFIG <- list(
    enabled = TRUE,           # TRUE for testing single comparison
    comparison = "comp_1108forNoneAndWT",  # Which comparison to process in debug mode
    save_plots = FALSE,        # Whether to save plots to files
    verbose = TRUE,           # Print debug information
    chromosome = 10,
    interactive = FALSE,
    validate_config = TRUE,
    display_time = 2
)

# Time formatting configuration
TIME_CONFIG <- list(
    timestamp_format = "%Y%m%d_%H%M%S",  # YYYYMMDD_HHMMSS
    date_format = "%Y%m%d"               # YYYYMMDD
)

# Generate timestamps once at script start
TIMESTAMPS <- list(
    full = format(Sys.time(), TIME_CONFIG$timestamp_format),
    date = format(Sys.Date(), TIME_CONFIG$date_format)
)

PLOT_CONFIG <- list(
    width = 10,
    height = 8,
    placeholder_color = "#cccccc",
    input_color ="#808080",
    track_name_format = "%s: %s",
    control_track_name_format = "%s: %s - %s",
    placeholder_suffix = "(No data)",

    #Title configuration
    title = list(
        mode = "development",  # or "publication"
        development = list(
            width = 0.9,
            fontface = 1,
            cex = 0.75,
            background = "white",
            # Format title with visual spacing
            format = paste(
                "%s",
                "Comparison: %s",
                "Chromosome %s (%d samples)",
                "%s",
                "Normalization: %s",
                sep = "\n"
            )
        ),
        publication = list(
            width = 0.4,
            fontface = 2,
            cex = 1,
            background = "transparent",
            format = "%s: Chr%s (%s)"
        ),
        format = list(
            max_width = 40,
            max_lines = 5
        )
    )

)

# Load required packages
#-----------------------------------------------------------------------------
required_packages <- c("rtracklayer", "GenomicRanges", "Gviz")
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' is missing", pkg))
    }
}

#source("~/lab_utils/failsafe_scripts/all_functions.R")
source("~/lab_utils/failsafe_scripts/functions_for_genome_tracks.R")
source("~/lab_utils/failsafe_scripts/bmc_config.R")
if (DEBUG_CONFIG$validate_config) {
    if (!exists("EXPERIMENT_CONFIG") ||
        !("COMPARISONS" %in% names(EXPERIMENT_CONFIG))) {
        stop("EXPERIMENT_CONFIG must contain COMPARISONS element")
    }

    if (DEBUG_CONFIG$verbose) {
        message("Found ", length(EXPERIMENT_CONFIG$COMPARISONS),
                " comparisons in configuration")
    }
}


# Load metadata and files
#-----------------------------------------------------------------------------
#experiment_id <- "241007Bel"
experiment_id <- "241010Bel"
base_dir <- file.path(Sys.getenv("HOME"), "data", experiment_id)
plots_dir <- file.path(base_dir, "plots", "genome_tracks", "comparisons")
metadata_path <- file.path(base_dir, "documentation",
                          paste0(experiment_id, "_sample_grid.csv"))
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

if (DEBUG_CONFIG$verbose) {
    message("Metadata processing summary:")
    message(sprintf("Found %d fastq files", length(fastq_files)))
    message(sprintf("Processed %d metadata rows", nrow(sorted_metadata)))
    message("Columns with enforced factors:")
    print(names(EXPERIMENT_CONFIG$CATEGORIES))
}

# Create color scheme
color_scheme <- create_color_scheme(
    config = list(
        placeholder = PLOT_CONFIG$placeholder_color,
        input = PLOT_CONFIG$input_color
    ),
    categories = list(
        antibody = unique(sorted_metadata$antibody),
        rescue_allele = unique(sorted_metadata$rescue_allele)
    )
)

if (DEBUG_CONFIG$verbose) {
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
chromosome_to_plot <- DEBUG_CONFIG$chromosome
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
pattern_for_feature_file <- "eaton_peaks"
feature_file <- list.files(
    file.path(Sys.getenv("HOME"), "data", "feature_files"),
    pattern = pattern_for_feature_file,
    full.names = TRUE
)[1]
# Convert to title case and replace "_" with " "
feature_track_name <- gsub(
                        "_",
                        " ",
                        tools::toTitleCase(pattern_for_feature_file)
                        )

if (!is.null(feature_file)) {
    features <- rtracklayer::import(feature_file)
    # Convert to chrRoman format
    GenomeInfoDb::seqlevels(features) <- paste0(
        "chr",
        utils::as.roman(gsub("chr", "", GenomeInfoDb::seqlevels(features)))
    )
}


# Process Comparisons
#-----------------------------------------------------------------------------
# Determine which comparisons to process
comparisons_to_process <- if (DEBUG_CONFIG$enabled) {
    if (!DEBUG_CONFIG$comparison %in% names(EXPERIMENT_CONFIG$COMPARISONS)) {
        stop("Debug comparison not found in EXPERIMENT_CONFIG$COMPARISONS")
    }
    DEBUG_CONFIG$comparison
} else {
    names(EXPERIMENT_CONFIG$COMPARISONS)
}

if (DEBUG_CONFIG$verbose) {
    message("\nProcessing comparisons:")
    message(sprintf("- Total comparisons: %d", length(comparisons_to_process)))
    message("- Comparisons: ", paste(comparisons_to_process, collapse = ", "))
}

# Calculate global y-limits for all plots (before the plotting loop)
#-----------------------------------------------------------------------------
if (DEBUG_CONFIG$verbose) {
    message("\nCalculating global range for all tracks...")
}

limits_result <- calculate_track_limits(
    bigwig_files = bigwig_files,
    genome_range = genome_range,
    padding_fraction = 0.1,
    verbose = DEBUG_CONFIG$verbose
)

if (!limits_result$success) {
    warning("Failed to calculate y-limits: ", limits_result$error)
    y_limits <- c(0, 1000)  # Default fallback
} else {
    y_limits <- limits_result$data
}

# Add after processing metadata but before track creation
short_sample_ids <- create_minimal_identifiers(
    sorted_metadata$sample_id,
    verbose = DEBUG_CONFIG$verbose
)

# Create mapping between full and short IDs
sample_id_mapping <- setNames(short_sample_ids, sorted_metadata$sample_id)

# Process each comparison
#-----------------------------------------------------------------------------
for (comparison_name in comparisons_to_process) {
    if (DEBUG_CONFIG$verbose) {
        message(sprintf("\nProcessing comparison: %s", comparison_name))
    }

    # Get comparison expression and filter metadata
    comparison_expression <- EXPERIMENT_CONFIG$COMPARISONS[[comparison_name]]
    comparison_samples <- sorted_metadata[eval(comparison_expression,
                                             envir = sorted_metadata), ]
    label_result <- create_track_labels(
        samples = comparison_samples,
        always_show = "antibody",
        categories = EXPERIMENT_CONFIG$CONTROL_FACTORS$genotype,
        verbose = DEBUG_CONFIG$verbose
    )

    if (!label_result$success) {
        warning("Failed to create track labels for comparison %s: %s", comparison_name, label_result$error)
        track_labels <- comparison_samples$short_name  # Fallback to sample_id
    } else {
        track_labels <- label_result$data$labels
    }

    if (nrow(comparison_samples) == 0) {
        warning(sprintf("No samples found for comparison: %s", comparison_name))
        next
    }

    if (DEBUG_CONFIG$verbose) {
        message(sprintf("Found %d samples for comparison", nrow(comparison_samples)))
    }

    # Initialize track list with genome axis
    tracks <- list(
        Gviz::GenomeAxisTrack(
            name = paste("Chr ", chromosome_to_plot, " Axis", sep = "")
        )
    )

    # Find and add control track (single control per comparison)
    if (DEBUG_CONFIG$verbose) {
        message("\nFinding control sample for comparison...")
    }

    # Try to find control for first available sample
    control_sample <- find_control_sample(
        comparison_samples[1, ],
        sorted_metadata,
        EXPERIMENT_CONFIG$CONTROL_FACTORS
    )

    # Add single control track (real or placeholder)
    if (!is.null(control_sample)) {
        control_bigwig <- bigwig_files[grepl(control_sample$sample_id,
                                           bigwig_files)]
        if (length(control_bigwig) > 0 && file.exists(control_bigwig[1])) {
            if (DEBUG_CONFIG$verbose) {
                message(sprintf("Adding control track: %s",
                              control_sample$sample_id))
            }

            control_track_data <- rtracklayer::import(
                control_bigwig[1],
                which = genome_range
            )

            tracks[[length(tracks) + 1]] <- Gviz::DataTrack(
                control_track_data,
                name = sprintf(
                    PLOT_CONFIG$control_track_name_format,
                    sample_id_mapping[control_sample$sample_id],
                    "Input",
                    control_sample$rescue_allele
                ),
                type = "l",
                col = color_scheme$fixed$input,
                chromosome = chromosome_roman,
                cex.title = 1.2
            )
        } else {
            if (DEBUG_CONFIG$verbose) {
                message("Creating placeholder for control track")
            }

            # Create placeholder control track
            empty_ranges <- GenomicRanges::GRanges(
                seqnames = chromosome_roman,
                ranges = IRanges::IRanges(
                    start = seq(1, chromosome_width, length.out = 1000),
                    width = 1
                ),
                score = rep(0, 1000),
                strand = "*"
            )

            tracks[[length(tracks) + 1]] <- Gviz::DataTrack(
                empty_ranges,
                name = sprintf(
                    "%s: %s - %s %s",
                    control_sample$sample_id,
                    "Input",
                    control_sample$rescue_allele,
                    PLOT_CONFIG$placeholder_suffix
                ),
                type = "l",
                col = color_scheme$fixed$placeholder,
                chromosome = chromosome_roman,
                cex.title = 1.2
            )
        }
    } else {
        if (DEBUG_CONFIG$verbose) {
            message("No matching control found, adding placeholder track")
        }

        # Create placeholder track when no control found
        empty_ranges <- GenomicRanges::GRanges(
            seqnames = chromosome_roman,
            ranges = IRanges::IRanges(
                start = seq(1, chromosome_width, length.out = 1000),
                width = 1
            ),
            score = rep(0, 1000),
            strand = "*"
        )

        tracks[[length(tracks) + 1]] <- Gviz::DataTrack(
            empty_ranges,
            name = sprintf("Input Control %s", PLOT_CONFIG$placeholder_suffix),
            type = "l",
            col = color_scheme$fixed$placeholder,
            chromosome = chromosome_roman,
            cex.title = 1.2
        )
    }

    # Add sample tracks
    for (i in seq_len(nrow(comparison_samples))) {
        sample_id <- comparison_samples$sample_id[i]
        sample_bigwig <- bigwig_files[grepl(sample_id, bigwig_files)]
        current_antibody <- comparison_samples$antibody[i]

        track_name <- sprintf(
            PLOT_CONFIG$track_name_format,
            sample_id_mapping[sample_id],
            track_labels[i]
        )

        if (length(sample_bigwig) > 0 && file.exists(sample_bigwig[1])) {
            if (DEBUG_CONFIG$verbose) {
                message(sprintf("Adding track for sample: %s", sample_id))
            }

            # Update track color based on antibody
            track_color <- if (current_antibody == "Input") {
                color_scheme$fixed$input
            } else {
                color_scheme$get_color("antibody", current_antibody)
            }

            track_data <- rtracklayer::import(
                sample_bigwig[1],
                which = genome_range
            )

            tracks[[length(tracks) + 1]] <- Gviz::DataTrack(
                track_data,
                name = track_name,
                type = "l",
                col = track_color,
                chromosome = chromosome_roman,
                cex.title = 1.2
            )
        } else {
            if (DEBUG_CONFIG$verbose) {
                message(sprintf("Creating placeholder for sample: %s", sample_id))
            }

            # Create placeholder track
            empty_ranges <- GenomicRanges::GRanges(
                seqnames = chromosome_roman,
                ranges = IRanges::IRanges(
                    start = seq(1, chromosome_width, length.out = 1000),
                    width = 1
                ),
                score = rep(0, 1000),
                strand = "*"
            )

            tracks[[length(tracks) + 1]] <- Gviz::DataTrack(
                empty_ranges,
                name = paste(track_name, PLOT_CONFIG$placeholder_suffix),
                type = "l",
                col = color_scheme$fixed$placeholder,
                chromosome = chromosome_roman,
                cex.title = 1.2
            )
        }
    }

    # Add feature track if available
    if (!is.null(features)) {
        tracks[[length(tracks) + 1]] <- Gviz::AnnotationTrack(
            features,
            name = feature_track_name,
            rotation.title = 0,
            cex.title = 1.2
        )
    }

    # For creating plot title, use category information
    distinguishing_categories <- label_result$data$categories$distinguishing
    varying_categories <- label_result$data$categories$varying
    all_categories <- label_result$data$categories$all_available


    # Find shared characteristics
    shared_categories <- setdiff(
        all_categories,
        c(distinguishing_categories, "sample_id")
    )


    # Create shared characteristics text
    shared_values <- sapply(shared_categories, function(col) {
        values <- unique(comparison_samples[[col]])
        if (length(values) == 1) {
            sprintf("%s: %s", col, values)
        } else {
            NULL
        }
    })
    shared_values <- unlist(shared_values[!sapply(shared_values, is.null)])

    # Debug output if needed
    if (DEBUG_CONFIG$verbose) {
        message("\nLabel Result Summary:")
        message("- Distinguishing categories: ",
                paste(distinguishing_categories, collapse = ", "))
        message("- Varying categories: ",
                paste(varying_categories, collapse = ", "))
        message("- Shared categories: ",
                paste(shared_categories, collapse = ", "))
        message("- Number of labels: ", length(track_labels))
    }

    # Create plot title
    #-----------------------------------------------------------------------------
    # Get title configuration for current mode
    title_config <- PLOT_CONFIG$title[[PLOT_CONFIG$title$mode]]

    plot_title <- sprintf(
        title_config$format,
        experiment_id,
        sub("^comp_", "", comparison_name),
        chromosome_to_plot,
        nrow(comparison_samples),
        format(Sys.time(), "%Y-%m-%d %H:%M"),
        normalization_method
    )

    # Generate plot
    #-----------------------------------------------------------------------------
    if (DEBUG_CONFIG$verbose) {
        message("\nGenerating visualization...")
    }


    Gviz::plotTracks(
        trackList = tracks,
        chromosome = chromosome_roman,
        from = 1,
        to = chromosome_width,
        ylim = y_limits,
        main = plot_title,
        # Title appearance based on mode
        title.width = title_config$width,
        fontface.title = title_config$fontface,
        cex.title = title_config$cex,
        background.title = title_config$background,
        col.title = "black",          # Title text color

        # Track name appearance
        fontcolor = "black",           # Track name text color
        background.title = "white",    # Track name background
        col.border.title = "#E0E0E0",  # Light gray border around track names
        fontsize = 11, # Base font size

        # Other visualization parameters
        cex.main = .75,
        fontface.main = 0.6,
        margin = 15,
        innerMargin = 5,

        # Axis appearance
        col.axis = "black", # Axis text color
        cex.axis = 0.8 # Axis text size

    )

    # Save plot if configured
    if (DEBUG_CONFIG$save_plots) {
        plot_file <- file.path(
            plots_dir,
            sprintf(
                "%s_%s_chr%s_%s.svg",
                TIMESTAMPS$full,
                experiment_id,
                chromosome_to_plot,
                sub("^comp_", "", comparison_name)
            )
        )

        if (DEBUG_CONFIG$verbose) {
            message(sprintf("Saving plot to: %s", basename(plot_file)))
        }

        svg(plot_file,
            width = PLOT_CONFIG$width,
            height = PLOT_CONFIG$height)
        Gviz::plotTracks(
            trackList = tracks,
            chromosome = chromosome_roman,
            from = 1,
            to = chromosome_width,
            ylim = y_limits,
            main = plot_title,
            # Title appearance based on mode
            title.width = title_config$width,
            rotation.title = title_config$rotation,
            fontface.title = title_config$fontface,
            cex.title = title_config$cex,
            background.title = title_config$background,

            # Title panel
            col.title = "black",          # Title text color
            fontface.title = 2,            # Bold title

            # Track name appearance
            fontcolor = "black",           # Track name text color
            background.title = "white",    # Track name background
            col.border.title = "#E0E0E0",  # Light gray border around track names

            # Other visualization parameters
            cex.main = 0.6,
            margin = 15,
            innerMargin = 5,

            # Axis appearance
            col.axis = "black",            # Axis text color
            cex.axis = 0.8               # Axis text size

        )
        dev.off()
    }

    # Interactive viewing options
    if (DEBUG_CONFIG$interactive) {
        user_input <- readline(
            prompt = "Options: [Enter] next comparison, 's' skip rest, 'q' quit: "
        )
        if (user_input == "q") break
        if (user_input == "s") DEBUG_CONFIG$save_plots <- FALSE
    } else {
        Sys.sleep(DEBUG_CONFIG$display_time)
    }
}


# Final processing summary
#-----------------------------------------------------------------------------
if (DEBUG_CONFIG$verbose) {
    message("\nProcessing complete")
    message(sprintf("Processed %d comparisons", length(comparisons_to_process)))
    if (DEBUG_CONFIG$save_plots) {
        message(sprintf("Output directory: %s", plots_dir))
    }
}
