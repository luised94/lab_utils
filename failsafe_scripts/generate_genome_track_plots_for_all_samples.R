#!/usr/bin/env Rscript
# Constants and Configuration
#-----------------------------------------------------------------------------
PLOT_CONFIG <- list(
    SAMPLES_PER_PAGE = 4,
    DEFAULT_CHROMOSOME = 10,
    TRACK_COLOR = "#fd0036",
    WIDTH = 10,
    HEIGHT = 8
)
EXPERIMENT_ID <- "241007Bel"
TIMESTAMP <- format(Sys.Date(), "%Y%m%d")
ref_genome_file <- list.files(
    file.path(Sys.getenv("HOME"), "data", "REFGENS"),
    pattern = "S288C_refgenome.fna",
    full.names = TRUE,
    recursive = TRUE
)[1]
feature_file <- list.files(
    file.path(Sys.getenv("HOME"), "data", "feature_files"),
    pattern = "eaton_peaks",
    full.names = TRUE
)[1]
REQUIRED_PACKAGES <- c("ShortRead", "GenomeInfoDb", "rtracklayer", "GenomicRanges", "Gviz", "tidyverse", "QuasR")

# 1. Load and verify config and functions
#-----------------------------------------------------------------------------
source("~/lab_utils/failsafe_scripts/all_functions.R")
source("~/lab_utils/failsafe_scripts/bmc_config.R")

#config_result <- validate_dependencies()
#if (!config_result) {
#    stop("Configuration validation failed")
#}
## 2. Load and verify required packages
#-----------------------------------------------------------------------------
packages_result <- packages_required_validate(REQUIRED_PACKAGES)
if (!packages_result$success) {
    stop(packages_result$error)
}

# 3. Set up initial variables and paths
#-----------------------------------------------------------------------------
experiment_result <- experiment_environment_validate(
    EXPERIMENT_ID,
    list(directories = c("coverage", "documentation", "plots"))
)
if (!experiment_result$success) {
    stop(experiment_result$error)
}

base_directory <- experiment_result$data$base_path
coverage_directory <- file.path(base_directory, "coverage")
plots_directory <- file.path(base_directory, "plots")


# 4. Load and process sample table
#-----------------------------------------------------------------------------
metadata_path_result <- metadata_path_validate(
    base_directory,
    "%s_sample_grid.csv"
)
if (!metadata_path_result$success) {
    stop(metadata_path_result$error)
}

metadata_result <- metadata_file_read(
    metadata_path_result$data,
    list(stringsAsFactors = FALSE)
)
if (!metadata_result$success) {
    stop(metadata_result$error)
}

# Enforce factor levels and sort
sample_table_processed <- factor_categories_enforce(
    metadata_result$data,
    EXPERIMENT_CONFIG$CATEGORIES
)$data

sample_table_sorted <- metadata_frame_sort(
    sample_table_processed,
    EXPERIMENT_CONFIG$COLUMN_ORDER
)$data

# 5. Load reference genome
#-----------------------------------------------------------------------------
genome_result <- tryCatch({
    Biostrings::readDNAStringSet(ref_genome_file)
}, error = function(e) {
    stop(sprintf("Failed to load reference genome: %s", e$message))
})

# 6. Create genome ranges
#-----------------------------------------------------------------------------
chromosome <- PLOT_CONFIG$DEFAULT_CHROMOSOME
genome_range_result <- genomic_range_create(
    chromosome,
    list(start = 1, end = 1e6)
)
if (!genome_range_result$success) {
    stop(genome_range_result$error)
}

# 7. Load and process feature file
#-----------------------------------------------------------------------------
feature_result <- tryCatch({
    features <- rtracklayer::import(feature_file)
    
    # Convert chromosome style
    GenomeInfoDb::seqlevels(features) <- chromosome_names_convert(
        GenomeInfoDb::seqlevels(features),
        "Roman"
    )$data
    
    feature_track <- feature_track_create(
        features,
        list(
            name = "Features",
            chromosome = genome_range_result$data@seqnames[1]
        )
    )
    
    if (!feature_track$success) {
        stop(feature_track$error)
    }
    
    feature_track$data
}, error = function(e) {
    warning(sprintf("Failed to load feature file: %s", e$message))
    NULL
})

# Process all samples
#-----------------------------------------------------------------------------
# Get all bigwig files
bigwig_files <- list.files(
    coverage_directory,
    pattern = ".*normalized.*\\.bw$",
    full.names = TRUE
)

# Calculate global range
range_result <- get_global_range(bigwig_files, genome_range_result$data)
if (!range_result$success) {
    stop(range_result$error)
}

# Process samples in groups
sample_groups <- split(
    seq_len(nrow(sample_table_sorted)),
    ceiling(seq_len(nrow(sample_table_sorted)) / PLOT_CONFIG$SAMPLES_PER_PAGE)
)

# Create plots for each group
for (group_idx in seq_along(sample_groups)) {
    group_samples <- sample_table_sorted[sample_groups[[group_idx]], ]
    
    # Create tracks for group
    track_group_result <- track_group_create(
        lapply(seq_len(nrow(group_samples)), function(i) {
            list(
                bigwig_file = bigwig_files[i],
                name = group_samples$sample_id[i]
            )
        }),
        list(
            chromosome = genome_range_result$data@seqnames[1],
            color = PLOT_CONFIG$TRACK_COLOR
        )
    )
    
    if (!track_group_result$success) {
        warning(sprintf("Failed to create tracks for group %d: %s", 
                       group_idx, track_group_result$error))
        next
    }
    
    # Add feature track if available
    if (!is.null(feature_result)) {
        track_group_result$data$tracks <- c(
            track_group_result$data$tracks,
            feature_result$track
        )
    }
    
    # Generate plot
    plot_path_result <- plot_path_generate(
        plots_directory,
        list(
            chromosome = chromosome,
            group = group_idx
        )
    )
    
    if (!plot_path_result$success) {
        warning(sprintf("Failed to generate plot path for group %d: %s",
                       group_idx, plot_path_result$error))
        next
    }
    
    plot_result <- plot_tracks_create(
        track_group_result$data,
        list(
            width = PLOT_CONFIG$WIDTH,
            height = PLOT_CONFIG$HEIGHT,
            output_path = plot_path_result$data,
            calculate_limits = TRUE,
            chromosome = genome_range_result$data@seqnames[1]
        )
    )
    
    if (!plot_result$success) {
        warning(sprintf("Failed to create plot for group %d: %s",
                       group_idx, plot_result$error))
        next
    }
    
    message(sprintf("Created plot for group %d: %s",
                   group_idx, basename(plot_result$data$path)))
}

message("Processing complete")
# 1. Load required packages with verification
# 2. Set up initial variables and paths
# 3. Load and verify config and functions
# 4. Load and process sample table with new processing steps
# Enforce factor levels from config
# Sort metadata using config column order
# 5. Load reference genome
# 6. Create genome ranges
# 7. Load feature file
# Load and adjust feature chromosome style
# Function for bigwig validation
# Setup chromosome
# Determine the maximum and minimum across all sample tracks. 
# Plot all genome tracks samples four at time using the previously determined maximum and minimum.
