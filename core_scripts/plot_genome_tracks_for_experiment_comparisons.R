#!/usr/bin/env Rscript
# replace all sorted_metadata
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
is_template <- args$is_template

file_directory <- if (is_template) args$experiment_dir else file.path(args$experiment_dir, "documentation")
file_identifier <- if (is_template) "template" else args$experiment_id

config_path <- file.path(file_directory, paste0(file_identifier, "_bmc_config.R"))
metadata_path <- file.path(file_directory, paste0(file_identifier, "_sample_grid.csv"))


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

if(!is.null(args$output_format)) {
    RUNTIME_CONFIG$output_format <- args$output_format

    message(sprintf("Setting output format: %s", RUNTIME_CONFIG$output_format))
}

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

reference_bigwig <- "~/data/100303Bel/coverage/processed_034475_sequence_to_S288C_CPM.bw"
if (!is.null(reference_bigwig)) {
    reference_grange <- rtracklayer::import(reference_bigwig)
    #print(reference_grange)
    ## Convert to chrRoman format
    #GenomeInfoDb::seqlevels(features) <- paste0(
    #    "chr",
    #    utils::as.roman(gsub("chr", "", GenomeInfoDb::seqlevels(features)))
    #)
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
        ".Target" = chromosome_to_plot,
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
comparisons_to_process <- if (RUNTIME_CONFIG$process_single_file) {
    if (!RUNTIME_CONFIG$process_comparison %in% names(EXPERIMENT_CONFIG$COMPARISONS)) {
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

for (comparison_name in comparisons_to_process) {
    # Metadata processing and sample selection
    comparison_expression <- EXPERIMENT_CONFIG$COMPARISONS[[comparison_name]]
    row_samples_to_visualize <- sorted_metadata[eval(comparison_expression,
                                             envir = sorted_metadata), ]

    # Try to find control for first available sample
    control_sample <- find_control_sample(
        row_samples_to_visualize[1, ],
        sorted_metadata,
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
           # name = paste("Chr ", chromosome_to_plot, " Axis", sep = "")
            name = sprintf(track_axis_name_template, chromosome_to_plot)
        )
    )

    # Try to find control for first available sample
    control_sample <- find_control_sample(
        row_samples_to_visualize[1, ],
        sorted_metadata,
        EXPERIMENT_CONFIG$CONTROL_FACTORS
    )
    
    if(!is.null(control_sample)) {
        sample_id <- control_sample$sample_id[i]
        current_antibody <- control_sample$antibody[i]
        track_label <- track_labels[i]
        bigwig_file_path <- project_bigwig_files[grepl(sample_id, project_bigwig_files)][1]
        control_bigwig_file_path <- bigwig_files[grepl(control_sample$sample_id,
                                           bigwig_files)]

    }
