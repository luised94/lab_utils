
#!/usr/bin/env Rscript

# Constants and Configuration
#-----------------------------------------------------------------------------
PLOT_CONFIG <- list(
    SAMPLES_PER_PAGE = 4,
    DEFAULT_CHROMOSOME = 10,
    TRACK_COLOR = "#fd0036",
    PLACEHOLDER_COLOR = "#cccccc",
    WIDTH = 10,
    HEIGHT = 8
)

EXPERIMENT_ID <- "241007Bel"
TIMESTAMP <- format(Sys.Date(), "%Y%m%d")

# Define file paths with intention-revealing names
reference_genome_path <- list.files(
    path = file.path(Sys.getenv("HOME"), "data", "REFGENS"),
    pattern = "S288C_refgenome.fna",
    full.names = TRUE,
    recursive = TRUE
)[1]

feature_annotation_path <- list.files(
    path = file.path(Sys.getenv("HOME"), "data", "feature_files"),
    pattern = "eaton_peaks",
    full.names = TRUE
)[1]

options(ucscChromosomeNames = FALSE)

REQUIRED_PACKAGES <- c(
    "ShortRead", "GenomeInfoDb", "rtracklayer", 
    "GenomicRanges", "Gviz", "tidyverse", "QuasR"
)

# 1. Load and verify dependencies
#-----------------------------------------------------------------------------
source("~/lab_utils/failsafe_scripts/all_functions.R")
source("~/lab_utils/failsafe_scripts/bmc_config.R")

packages_validation_result <- packages_required_validate(REQUIRED_PACKAGES)
if (!packages_validation_result$success) {
    stop(packages_validation_result$error)
}

# 2. Initialize experiment environment
#-----------------------------------------------------------------------------
experiment_validation_result <- experiment_environment_validate(
    experiment_identifier = EXPERIMENT_ID,
    environment_requirements = list(
        directories = c("coverage", "documentation", "plots")
    )
)
if (!experiment_validation_result$success) {
    stop(experiment_validation_result$error)
}

experiment_paths <- list(
    base = experiment_validation_result$data$base_path,
    coverage = file.path(experiment_validation_result$data$base_path, "coverage"),
    plots = file.path(experiment_validation_result$data$base_path, "plots")
)

# 3. Load and process metadata
#-----------------------------------------------------------------------------
metadata_path_result <- metadata_path_validate(
    directory_path = experiment_paths$base,
    file_pattern = "%s_sample_grid.csv"
)
if (!metadata_path_result$success) {
    stop(metadata_path_result$error)
}

metadata_result <- metadata_file_read(
    file_path = metadata_path_result$data,
    read_options = list(stringsAsFactors = FALSE)
)
if (!metadata_result$success) {
    stop(metadata_result$error)
}

# Process metadata with factor levels and sorting
processed_metadata <- factor_categories_enforce(
    metadata_frame = metadata_result$data,
    category_definitions = EXPERIMENT_CONFIG$CATEGORIES
)$data

sorted_metadata <- metadata_frame_sort(
    metadata_frame = processed_metadata,
    sort_columns = EXPERIMENT_CONFIG$COLUMN_ORDER
)$data

# 4. Load genome and feature data
#-----------------------------------------------------------------------------
genome_data <- tryCatch({
    Biostrings::readDNAStringSet(reference_genome_path)
}, error = function(e) {
    stop(sprintf("Failed to load reference genome: %s", e$message))
})

genome_range_result <- genomic_range_create(
    chromosome_number = PLOT_CONFIG$DEFAULT_CHROMOSOME,
    range_parameters = list(start = 1, end = 1e6)
)
if (!genome_range_result$success) {
    stop(genome_range_result$error)
}

# Load and process feature annotations
feature_track_result <- tryCatch({
    feature_data <- rtracklayer::import(feature_annotation_path)
    
    # Convert chromosome style
    GenomeInfoDb::seqlevels(feature_data) <- chromosome_names_convert(
        chromosome_names = GenomeInfoDb::seqlevels(feature_data),
        target_style = "Roman"
    )$data
    
    feature_track <- feature_track_create(
        feature_data = feature_data,
        track_options = list(
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

# 5. Process bigwig files and create sample mapping
#-----------------------------------------------------------------------------
available_bigwig_files <- list.files(
    path = experiment_paths$coverage,
    pattern = "_CPM\\.bw$",
    full.names = TRUE
)

bigwig_mapping_result <- create_bigwig_sample_mapping(
    sample_table = sorted_metadata,
    bigwig_files = available_bigwig_files
)
if (!bigwig_mapping_result$success) {
    stop(bigwig_mapping_result$error)
}

# Calculate global range for visualization
range_result <- get_global_range(
    bigwig_files = available_bigwig_files,
    genome_range = genome_range_result$data
)
if (!range_result$success) {
    stop(range_result$error)
}

# 6. Process samples in groups and create plots
#-----------------------------------------------------------------------------
# Calculate grouping parameters
sample_count <- nrow(sorted_metadata)
group_count <- ceiling(sample_count / PLOT_CONFIG$SAMPLES_PER_PAGE)

# Create sample groups
sample_groups <- split(
    x = seq_len(sample_count),
    f = ceiling(seq_len(sample_count) / PLOT_CONFIG$SAMPLES_PER_PAGE)
)

# Process each group
for (group_idx in seq_along(sample_groups)) {
    current_group_samples <- sorted_metadata[sample_groups[[group_idx]], ]
    
    # Create track configurations for current group
    track_configs_result <- create_sample_track_configs(
        group_samples = current_group_samples,
        bigwig_mapping = bigwig_mapping_result$data
    )
    if (!track_configs_result$success) {
        warning(sprintf("Failed to create track configs for group %d: %s",
                       group_idx, track_configs_result$error))
        next
    }
    
    # Create visualization tracks
    track_group_result <- track_group_create(
        sample_list = track_configs_result$data,
        group_options = list(
            chromosome = genome_range_result$data@seqnames[1],
            color = PLOT_CONFIG$TRACK_COLOR,
            placeholder_color = PLOT_CONFIG$PLACEHOLDER_COLOR
        )
    )
    if (!track_group_result$success) {
        warning(sprintf("Failed to create tracks for group %d: %s",
                       group_idx, track_group_result$error))
        next
    }
    
    # Add feature track if available
    if (!is.null(feature_track_result)) {
        track_group_result$data$tracks <- c(
            track_group_result$data$tracks,
            feature_track_result$track
        )
    }
    
    # Generate plot path
    plot_path_result <- plot_path_generate(
        base_directory = experiment_paths$plots,
        plot_parameters = list(
            chromosome = PLOT_CONFIG$DEFAULT_CHROMOSOME,
            group = group_idx,
            timestamp = TIMESTAMP,
            experiment_id = EXPERIMENT_ID
        )
    )
    if (!plot_path_result$success) {
        warning(sprintf("Failed to generate plot path for group %d: %s",
                       group_idx, plot_path_result$error))
        next
    }
    
    # Create visualization plot
    plot_result <- plot_tracks_create(
        track_group = track_group_result$data,
        plot_options = list(
            width = PLOT_CONFIG$WIDTH,
            height = PLOT_CONFIG$HEIGHT,
            output_path = plot_path_result$data,
            calculate_limits = TRUE,
            chromosome = genome_range_result$data@seqnames[1],
            title = sprintf("Chromosome %s - Group %d", 
                          PLOT_CONFIG$DEFAULT_CHROMOSOME, group_idx)
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
