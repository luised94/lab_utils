#!/usr/bin/env Rscript

# Configuration
#-----------------------------------------------------------------------------
DEBUG_CONFIG <- list(
    enabled = TRUE,           # TRUE for testing single group, FALSE for all
    group = 1,               # Which group to process when in debug mode
    samples_per_group = 4,    # Samples per plot
    save_plots = FALSE,       # Whether to save plots to files
    verbose = TRUE           # Print debug information
)

PLOT_CONFIG <- list(
    width = 10,
    height = 8,
    track_color = "#fd0036",
    placeholder_color = "#cccccc"
)

# Load required packages
#-----------------------------------------------------------------------------
required_packages <- c("rtracklayer", "GenomicRanges", "Gviz")
for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' is missing", pkg))
    }
}

source("~/lab_utils/failsafe_scripts/all_functions.R")
source("~/lab_utils/failsafe_scripts/bmc_config.R")
# Load metadata and files
#-----------------------------------------------------------------------------
experiment_id <- "241007Bel"
base_dir <- file.path(Sys.getenv("HOME"), "data", experiment_id)
plots_dir <- file.path(base_dir, "plots", "genome_tracks")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# Load sample metadata
metadata_processing_result <- experiment_metadata_process(
    directory_path = base_dir,
    configuration = list(
        categories = EXPERIMENT_CONFIG$CATEGORIES,
        column_order = EXPERIMENT_CONFIG$COLUMN_ORDER
    ),
    output_options = list(
        output_file = FALSE,  # Set to TRUE if you want to save processed metadata
        output_path = NULL
    )
)
if (!metadata_processing_result$success) {
    stop(metadata_processing_result$error)
}
sorted_metadata <- metadata_processing_result$data
# Find bigwig files
bigwig_files <- list.files(
    file.path(base_dir, "coverage"),
    pattern = "_CPM\\.bw$",
    full.names = TRUE
)

# Load feature file (annotation)
feature_file <- list.files(
    file.path(Sys.getenv("HOME"), "data", "feature_files"),
    pattern = "eaton_peaks",
    full.names = TRUE
)[1]

if (!is.null(feature_file)) {
    features <- rtracklayer::import(feature_file)
    # Make sure chromosome naming matches bigwig files
    GenomeInfoDb::seqlevels(features) <- paste0("chr", GenomeInfoDb::seqlevels(features))
}

# Process samples in groups
#-----------------------------------------------------------------------------
# Create sample groups
sample_groups <- split(
    seq_len(nrow(sorted_metadata)),
    ceiling(seq_len(nrow(sorted_metadata)) / DEBUG_CONFIG$samples_per_group)
)

# Determine which groups to process
groups_to_process <- if (DEBUG_CONFIG$enabled) {
    DEBUG_CONFIG$group
} else {
    seq_along(sample_groups)
}

# Process each group
for (group_idx in groups_to_process) {
    if (DEBUG_CONFIG$verbose) {
        message("\nProcessing group ", group_idx)
    }
    
    # Get current group's samples
    current_samples <- sorted_metadata[sample_groups[[group_idx]], ]
    
    # Initialize tracks list with chromosome axis
    tracks <- list(Gviz::GenomeAxisTrack())
    
    # Add tracks for each sample
    for (i in seq_len(nrow(current_samples))) {
        sample_id <- current_samples$sample_id[i]
        
        # Find matching bigwig file
        bigwig_file <- bigwig_files[grepl(sample_id, bigwig_files)][1]
        
        if (!is.na(bigwig_file) && file.exists(bigwig_file)) {
            if (DEBUG_CONFIG$verbose) {
                message("Adding track for sample: ", sample_id)
            }
            
            track_data <- rtracklayer::import(bigwig_file)
            
            tracks[[length(tracks) + 1]] <- Gviz::DataTrack(
                track_data,
                name = sample_id,
                type = "l",
                col = PLOT_CONFIG$track_color
            )
        } else {
            if (DEBUG_CONFIG$verbose) {
                message("Creating placeholder for sample: ", sample_id)
            }
            
            # Simple placeholder track
            empty_track <- Gviz::DataTrack(
                GenomicRanges::GRanges(),
                name = paste(sample_id, "(No Data)"),
                type = "l",
                col = PLOT_CONFIG$placeholder_color
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
    
    # Generate plot
    if (DEBUG_CONFIG$save_plots) {
        plot_file <- file.path(
            plots_dir,
            sprintf("%s_group%d.svg", experiment_id, group_idx)
        )
        
        if (DEBUG_CONFIG$verbose) {
            message("Saving plot to: ", plot_file)
        }
        
        svg(plot_file, width = PLOT_CONFIG$width, height = PLOT_CONFIG$height)
        Gviz::plotTracks(tracks)
        dev.off()
    } else {
        Gviz::plotTracks(tracks)
    }
}

message("Processing complete")
