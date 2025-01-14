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
short_sample_ids <- create_minimal_identifiers(metadata$sample_id, verbose = RUNTIME_CONFIG$debug_verbose)
# Create mapping between full and short IDs
sample_id_mapping <- setNames(short_sample_ids, metadata$sample_id)

if (RUNTIME_CONFIG$debug_verbose) {
    message("\nColor scheme initialized:")
    message("Fixed colors:")
    print(color_scheme$fixed)
    message("Category colors:")
    print(color_scheme$categories)
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
