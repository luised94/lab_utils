#!/usr/bin/env Rscript

# Script: genome_track_visualization.R
# Purpose: Generate genome browser tracks from BigWig files with feature annotations
# Usage: Rscript genome_track_visualization.R
# Dependencies: rtracklayer, GenomicRanges, Gviz
# Author: [AUTHOR]
# Date: [DATE]

# Required Packages
#-----------------------------------------------------------------------------
required_packages <- c("rtracklayer", "GenomicRanges", "Gviz", "Biostrings")
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Required package '%s' is missing", pkg))
    }
}

# Load Custom Functions and Configurations
#-----------------------------------------------------------------------------
# MODIFY PATHS AS NEEDED
source("path/to/config.R")     # Load experiment configurations
#source("path/to/functions.R")  # Load custom functions if needed

# File Paths and Directory Setup
#-----------------------------------------------------------------------------
# MODIFY THESE VARIABLES FOR YOUR EXPERIMENT
experiment_id <- "EXPERIMENT_ID"  # Replace with your experiment ID

# Construct directory paths
paths <- list(
    base = file.path(Sys.getenv("HOME"), "data", experiment_id),
    plots = NULL,
    coverage = NULL,
    documentation = NULL,
    reference = file.path(Sys.getenv("HOME"), "data", "REFGENS"),
    features = file.path(Sys.getenv("HOME"), "data", "feature_files")
)

# Initialize derived paths
paths$plots <- file.path(paths$base, "plots", "genome_tracks")
paths$coverage <- file.path(paths$base, "coverage")
paths$documentation <- file.path(paths$base, "documentation")

# Create output directory
dir.create(paths$plots, recursive = TRUE, showWarnings = FALSE)

# Data Loading and Preprocessing
#-----------------------------------------------------------------------------
# Load and validate metadata
metadata_path <- file.path(
    paths$documentation, 
    sprintf("%s_sample_grid.csv", experiment_id)
)

if (!file.exists(metadata_path)) {
    stop("Metadata file not found: ", metadata_path)
}

# Load sample metadata
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)

# Extract sample IDs from fastq files
#-----------------------------------------------------------------------------
fastq_files <- list.files(
    path = file.path(paths$base, "fastq"),
    pattern = "consolidated_.*_sequence\\.fastq$",
    full.names = FALSE
)

if (length(fastq_files) == 0) {
    stop("No fastq files found in: ", file.path(paths$base, "fastq"))
}

# Extract sample IDs
sample_ids <- gsub(
    pattern = "consolidated_([0-9]{5,6})_sequence\\.fastq",
    replacement = "\\1",
    x = fastq_files
)

# Process Metadata
#-----------------------------------------------------------------------------
# Enforce factor levels from configuration
for (col_name in names(EXPERIMENT_CONFIG$CATEGORIES)) {
    if (col_name %in% colnames(metadata)) {
        metadata[[col_name]] <- factor(
            metadata[[col_name]],
            levels = EXPERIMENT_CONFIG$CATEGORIES[[col_name]],
            ordered = TRUE
        )
    }
}

# Sort metadata and add sample IDs
sorted_metadata <- metadata[do.call(order, metadata[EXPERIMENT_CONFIG$COLUMN_ORDER]), ]
sorted_metadata$sample_id <- sample_ids

# Print processing summary if verbose
if (RUNTIME_CONFIG$debug_verbose) {
    message("\nMetadata Processing Summary:")
    message(sprintf("- Found %d fastq files", length(fastq_files)))
    message(sprintf("- Processed %d metadata rows", nrow(sorted_metadata)))
    message("- Enforced factors for columns:")
    print(names(EXPERIMENT_CONFIG$CATEGORIES))
}

# Load Reference Genome and Create Range
#-----------------------------------------------------------------------------
# Find reference genome file
ref_genome_file <- list.files(
    paths$reference,
    pattern = "S288C_refgenome.fna",
    full.names = TRUE,
    recursive = TRUE
)[1]

if (is.na(ref_genome_file)) {
    stop("Reference genome file not found in: ", paths$reference)
}

# Load genome and create range
genome_data <- Biostrings::readDNAStringSet(ref_genome_file)
chromosome_to_plot <- RUNTIME_CONFIG$process_chromosome
chromosome_width <- genome_data[chromosome_to_plot]@ranges@width
chromosome_roman <- paste0("chr", utils::as.roman(chromosome_to_plot))

genome_range <- GenomicRanges::GRanges(
    seqnames = chromosome_roman,
    ranges = IRanges::IRanges(start = 1, end = chromosome_width),
    strand = "*"
)

# Load BigWig Files and Features
#-----------------------------------------------------------------------------
# Find bigwig files
bigwig_files <- list.files(
    paths$coverage,
    pattern = "_CPM\\.bw$",
    full.names = TRUE
)

if (length(bigwig_files) == 0) {
    warning("No BigWig files found in: ", paths$coverage)
}

# Load feature annotations if available
feature_file <- list.files(
    paths$features,
    pattern = "eaton_peaks",
    full.names = TRUE
)[1]

features <- NULL
if (!is.null(feature_file)) {
    features <- tryCatch({
        feature_data <- rtracklayer::import(feature_file)
        # Convert to chrRoman format
        GenomeInfoDb::seqlevels(feature_data) <- paste0(
            "chr", 
            utils::as.roman(gsub("chr", "", GenomeInfoDb::seqlevels(feature_data)))
        )
        feature_data
    }, error = function(e) {
        warning("Failed to load feature file: ", e$message)
        NULL
    })
}

# Sample Grouping and Processing Setup
#-----------------------------------------------------------------------------
# Create sample groups for visualization
sample_groups <- split(
    seq_len(nrow(sorted_metadata)),
    ceiling(seq_len(nrow(sorted_metadata)) / RUNTIME_CONFIG$process_samples_per_group)
)

# Determine processing scope based on debug configuration
    RUNTIME_CONFIG$process_group
} else {
    seq_along(sample_groups)
}

if (RUNTIME_CONFIG$debug_verbose) {
    message("\nProcessing Configuration:")
    message(sprintf("- Total groups: %d", length(sample_groups)))
    message(sprintf("- Processing groups: %s", 
                   paste(groups_to_process, collapse = ", ")))
}

# Calculate Global Y-axis Limits
#-----------------------------------------------------------------------------
if (RUNTIME_CONFIG$debug_verbose) {
    message("\nCalculating global range for tracks...")
}

# Collect all track values for scaling
all_track_values <- c()

# Process each bigwig file
for (bigwig_file in bigwig_files) {
    if (file.exists(bigwig_file)) {
        tryCatch({
            track_data <- rtracklayer::import(
                bigwig_file,
                which = genome_range
            )
            
            if (length(track_data) > 0) {
                values <- GenomicRanges::values(track_data)$score
                if (length(values) > 0) {
                    all_track_values <- c(all_track_values, values)
                }
            }
        }, error = function(e) {
            if (RUNTIME_CONFIG$debug_verbose) {
                message("Skipping ", basename(bigwig_file), ": ", e$message)
            }
        })
    }
}

# Calculate visualization limits with padding
y_limits <- if (length(all_track_values) > 0) {
    y_min <- min(all_track_values, na.rm = TRUE)
    y_max <- max(all_track_values, na.rm = TRUE)
    y_range <- y_max - y_min
    c(
        y_min - (y_range * 0.1),  # Add 10% padding
        y_max + (y_range * 0.1)
    )
} else {
    c(0, 1)  # Default limits if no data
}

if (RUNTIME_CONFIG$debug_verbose) {
    message(sprintf("Global y-limits: [%.2f, %.2f]", y_limits[1], y_limits[2]))
}

# Track Creation and Visualization
#-----------------------------------------------------------------------------
for (group_idx in groups_to_process) {
    if (RUNTIME_CONFIG$debug_verbose) {
        message(sprintf("\nProcessing group %d of %d", 
                       group_idx, max(groups_to_process)))
    }
    
    # Get current group samples
    current_samples <- sorted_metadata[sample_groups[[group_idx]], ]
    
    # Initialize track list with genome axis
    tracks <- list(
        Gviz::GenomeAxisTrack(
            name = paste("Chr ", chromosome_to_plot, " Axis", sep = "")
        )
    )
    
    # Process each sample in group
    for (i in seq_len(nrow(current_samples))) {
        sample_id <- current_samples$sample_id[i]
        track_name <- sprintf(
            GENOME_TRACK_CONFIG$format_track,
            sample_id,
            current_samples$short_name[i]
        )
        
        # Find matching bigwig file
        bigwig_file <- bigwig_files[grepl(sample_id, bigwig_files)][1]
        
        if (!is.na(bigwig_file) && file.exists(bigwig_file)) {
            if (RUNTIME_CONFIG$debug_verbose) {
                message("Creating track for: ", track_name)
            }
            
            track_data <- rtracklayer::import(
                bigwig_file,
                which = genome_range
            )
            
            tracks[[length(tracks) + 1]] <- Gviz::DataTrack(
                track_data,
                name = track_name,
                type = "l",
                col = GENOME_TRACK_CONFIG$track_color
            )
        } else {
            if (RUNTIME_CONFIG$debug_verbose) {
                message("Creating placeholder for: ", track_name)
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
                name = paste(track_name, GENOME_TRACK_CONFIG$format_placeholder),
                type = "l",
                col = GENOME_TRACK_CONFIG$placeholder_color,
                chromosome = chromosome_roman
            )
        }
    }
    
    # Add feature track if available
    if (!is.null(features)) {
        tracks[[length(tracks) + 1]] <- Gviz::AnnotationTrack(
            features,
            name = "Features"
        )
    }
    
    # Generate plot title
    plot_title <- sprintf(
        "Chromosome %s - Group %d of %d",
        chromosome_to_plot,
        group_idx,
        max(groups_to_process)
    )
    
    # Create visualization
    Gviz::plotTracks(
        tracks,
        chromosome = chromosome_roman,
        from = 1,
        to = chromosome_width,
        ylim = y_limits,
        title = plot_title
    )
    
    # Save plot if configured
    if (RUNTIME_CONFIG$output_save_plots) {
        plot_file <- file.path(
            paths$plots,
            sprintf(
                "%s_%s_chr%s_group%d.svg",
                TIMESTAMPS$full,
                experiment_id,
                chromosome_to_plot,
                group_idx
            )
        )
        
        if (RUNTIME_CONFIG$debug_verbose) {
            message(sprintf("Saving plot to: %s", basename(plot_file)))
        }
        
        svg(plot_file, width = GENOME_TRACK_CONFIG$width, height = GENOME_TRACK_CONFIG$display_height)
        Gviz::plotTracks(
            tracks,
            chromosome = chromosome_roman,
            from = 1,
            to = chromosome_width,
            ylim = y_limits,
            title = plot_title
        )
        dev.off()
    }
    
    # Handle interactive viewing
    if (RUNTIME_CONFIG$debug_interactive) {
        user_input <- readline(
            prompt = "Options: [Enter] next plot, 's' skip rest, 'q' quit: "
        )
        if (user_input == "q") break
        if (user_input == "s") RUNTIME_CONFIG$output_save_plots <- FALSE
    } else {
        Sys.sleep(RUNTIME_CONFIG$output_display_time)
    }
}

message("\nProcessing complete")
