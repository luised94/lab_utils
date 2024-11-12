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

# Debug and execution configuration
PROCESSING_CONFIG <- list(
    debug_mode = TRUE,  # Toggle between debug (single group) and full processing
    debug_group = 15,    # Which group to process in debug mode
    verbose = TRUE      # Enable detailed logging
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

#options(ucscChromosomeNames = FALSE)

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

metadata_processing_result <- experiment_metadata_process(
    directory_path = experiment_paths$base,
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

# Verify sample_id column exists
if (!"sample_id" %in% colnames(sorted_metadata)) {
    stop("Processed metadata missing required 'sample_id' column")
}

# 4. Load genome and feature data
#-----------------------------------------------------------------------------
genome_data <- tryCatch({
    Biostrings::readDNAStringSet(reference_genome_path)
}, error = function(e) {
    stop(sprintf("Failed to load reference genome: %s", e$message))
})
chromosome_length <- genome_data[PLOT_CONFIG$DEFAULT_CHROMOSOME]@ranges@width

genome_range_result <- genomic_range_create(
    chromosome_number = PLOT_CONFIG$DEFAULT_CHROMOSOME,
    range_parameters = list(start = 1, end = chromosome_length)
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

    # Convert chromosome style
    GenomeInfoDb::seqlevels(feature_data) <- chromosome_names_convert(
        chromosome_names = GenomeInfoDb::seqlevels(feature_data),
        target_style = "UCSC"
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

# Determine groups to process based on configuration
groups_to_process <- if (PROCESSING_CONFIG$debug_mode) {
    PROCESSING_CONFIG$debug_group
} else {
    seq_along(sample_groups)
}

# Validation and logging of processing configuration
message(sprintf("Processing mode: %s", 
                if(PROCESSING_CONFIG$debug_mode) "DEBUG (single group)" else "FULL"))
message(sprintf("Total groups: %d", length(sample_groups)))
message(sprintf("Groups to process: %s", 
                paste(groups_to_process, collapse = ", ")))


# Bigwig file validation
if (PROCESSING_CONFIG$verbose) {
    message("\nValidating bigwig files:")
    bigwig_validation <- lapply(available_bigwig_files, function(file) {
        file_size <- file.size(file)
        message(sprintf(
            "File: %s\n  Exists: %s\n  Size: %d bytes",
            basename(file),
            file.exists(file),
            file_size
        ))
        return(file_size > 0)
    })
    
    valid_bigwig_count <- sum(unlist(bigwig_validation))
    message(sprintf("\nValid bigwig files: %d/%d\n",
                   valid_bigwig_count, length(available_bigwig_files)))
}

# Must see if feature_track_result can be accessed correclty. I also could add import which argument to pass genome range.
# Process selected groups

# Process sample groups
#-----------------------------------------------------------------------------
for (group_idx in groups_to_process) {
    message(sprintf("\nProcessing group %d", group_idx))
    
    # Get current group samples
    current_group_samples <- sorted_metadata[sample_groups[[group_idx]], ]
    
    # Create track configurations
    track_configs_result <- create_sample_track_configs(
        group_samples = current_group_samples,
        bigwig_mapping = bigwig_mapping_result$data
    )
    
    if (!track_configs_result$success) {
        warning(sprintf("Failed to create track configs for group %d: %s",
                       group_idx, track_configs_result$error))
        next
    }
    
    # Debug output for configurations
    if (PROCESSING_CONFIG$verbose) {
        message("\nTrack configurations for group ", group_idx, ":")
        print(data.frame(
            sample_id = sapply(track_configs_result$data, `[[`, "name"),
            bigwig = sapply(track_configs_result$data, function(x) basename(x$bigwig_file)),
            stringsAsFactors = FALSE
        ))
    }
    
    # Initialize track list with genome axis
    tracks <- list(
        Gviz::GenomeAxisTrack(
            name = sprintf("%s", genome_range_result$data@seqnames[1])
        )
    )
    
    # Process each sample in group
    for (track_config in track_configs_result$data) {
        if (PROCESSING_CONFIG$verbose) {
            message("\nProcessing sample: ", track_config$name)
        }
        
        # Create track (real or placeholder)
        if (is.na(track_config$bigwig_file) || !file.exists(track_config$bigwig_file)) {
            if (PROCESSING_CONFIG$verbose) {
                message("Creating placeholder for: ", track_config$name)
            }
            
            # Create placeholder track with proper genome range
            empty_ranges <- GenomicRanges::GRanges(
                seqnames = genome_range_result$data@seqnames[1],
                ranges = IRanges::IRanges(
                    start = 1,
                    end = genome_data[PLOT_CONFIG$DEFAULT_CHROMOSOME]@ranges@width
                ),
            )
            
            current_track <- Gviz::DataTrack(
                empty_ranges,
                type = "l",
                name = paste(track_config$name, "(No Data)"),
                col = PLOT_CONFIG$PLACEHOLDER_COLOR,
                chromosome = genome_range_result$data@seqnames[1]
            )
        } else {
            if (PROCESSING_CONFIG$verbose) {
                message("Creating track from: ", basename(track_config$bigwig_file))
            }
            
            # Import data for specific chromosome range
            track_data <- rtracklayer::import(
                track_config$bigwig_file,
                which = genome_range_result$data
            )
            
            current_track <- Gviz::DataTrack(
                track_data,
                type = "l",
                name = track_config$name,
                col = PLOT_CONFIG$TRACK_COLOR,
                chromosome = genome_range_result$data@seqnames[1]
            )
        }
        
        tracks[[length(tracks) + 1]] <- current_track
    }
    
    # Add feature track if available
    if (!is.null(feature_track_result)) {
        tracks[[length(tracks) + 1]] <- feature_track_result$track
    }
    
    # Create output directory if needed
    output_directory <- file.path(experiment_paths$plots, "genome_tracks/overview")
    dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
    
    # Generate plot filename
    plot_filename <- sprintf(
        "%s_%s_chr%s_group%d.svg",
        TIMESTAMP,
        EXPERIMENT_ID,
        PLOT_CONFIG$DEFAULT_CHROMOSOME,
        group_idx
    )
    output_path <- file.path(output_directory, plot_filename)
    
    if (PROCESSING_CONFIG$verbose) {
        message("\nGenerating plot: ", basename(output_path))
    }
    
    # Create plot
    tryCatch({
        #svg(output_path, width = PLOT_CONFIG$WIDTH, height = PLOT_CONFIG$HEIGHT)
        Gviz::plotTracks(
            tracks,
            chromosome = genome_range_result$data@seqnames[1],
            title = sprintf("Chromosome %s - Group %d", 
                          PLOT_CONFIG$DEFAULT_CHROMOSOME, group_idx),
            ylim = range_result$data,
            from = 1,
            to = genome_data[PLOT_CONFIG$DEFAULT_CHROMOSOME]@ranges@width
        )
        #dev.off()
        message(sprintf("Successfully created plot: %s", basename(output_path)))
    }, error = function(e) {
        warning(sprintf("Failed to create plot for group %d: %s", group_idx, e$message))
    })
}
message("Processing complete")
