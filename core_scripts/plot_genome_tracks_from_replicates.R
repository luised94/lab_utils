#!/usr/bin/env Rscript
source(file.path(Sys.getenv("HOME"), "lab_utils", "core_scripts", "functions_for_logging.R"))
source(file.path(Sys.getenv("HOME"), "lab_utils", "core_scripts", "functions_for_script_control.R"))

################################################################################
# Handle script arguments
################################################################################
# Parse arguments and validate configurations
description <- "Plot tracks from different experiments in same plot."
args <- parse_common_arguments(description = description)

# Proceed if packages are installed. Can be disable.
required_packages <- c("rtracklayer", "GenomicRanges", "Gviz")
check_required_packages(required_packages, verbose = TRUE, skip_validation = args$skip_validation)

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

args_info <- list(
    title = "Script Configuration",
    "script.name" = get_script_name(),
    "script.description" = description
)
print_debug_info(modifyList(args_info, args))

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
#todo:Create function that sources but reassigns EXPERIMENT_CONFIG to another name. Would need to update all references or find some way to update that. However, the settings should remain the same for the most part and I can include a safe_source call inside the for loop just in case.
sapply(config_path[1], safe_source)

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


# Handle configuration override (independent)
if (!is.null(args$override)) {
    structured_log_info("Starting override")
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
# Setup logging if requested (independent)
#if (args$log_to_file) {
#    log_file <- setup_logging("chip_peak_analysis")
#    flog.appender(appender.file(log_file))
#    flog.threshold(INFO)
#}

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

#####################
# Load the metadata dataframe with all experiments.
#####################
project_id <- c()
for (path in config_path) {
    safe_source(path, verbose = FALSE)
    experiment_dir <- dirname(dirname(path))
    metadata_path <- file.path(
        experiment_dir,
        "documentation",
        paste0(EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID, "_sample_grid.csv")
    )
    required_directories <- c("fastq", "documentation", "coverage")
    dirs <- setup_experiment_dirs(experiment_dir = experiment_dir,
        output_dir_name = "plots",
        required_input_dirs = required_directories
    )
    project_id <- c(project_id, EXPERIMENT_CONFIG$METADATA$PROJECT_ID)
    print(length(unique(project_id)) == 1)
    output_dir <- file.path(Sys.getenv("HOME"), "data", unique(project_id), "plots", "genome_tracks", "replicates")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
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

    if (length(bigwig_files) == 0) {
        stop("No bigwig files found in specified directory")
    }
    normalization_method <- sub(".*_([^_]+)\\.bw$", "\\1",
        basename(bigwig_files[1]))

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
        debug_info <- list(
            "title" = "File Processing Debug Information",
            "Directories" = NULL,
            ".FASTQ" = dirs$fastq,
            ".Coverage" = dirs$coverage,
            "Pattern Matching" = NULL,
            ".FASTQ Pattern" = GENOME_TRACK_CONFIG$file_pattern,
            ".Bigwig Pattern" = bigwig_pattern,
            "File Discovery" = NULL,
            ".FASTQ Files" = sprintf("Found: %d", length(fastq_files))
        )

        # Add first few FASTQ files if any exist
        if (length(fastq_files) > 0) {
            for (i in seq_along(head(fastq_files, 3))) {
                debug_info[[sprintf("..FASTQ File %d", i)]] <- fastq_files[i]
            }
        }

        # Add Bigwig information
        debug_info[[".Bigwig Files"]] <- sprintf("Found: %d", length(bigwig_files))
        if (length(bigwig_files) > 0) {
            for (i in seq_along(head(bigwig_files, 3))) {
                debug_info[[sprintf("..Bigwig File %d", i)]] <- basename(bigwig_files[i])
            }
            debug_info[[".Normalization"]] <- normalization_method
        }

        # Add Sample IDs information
        debug_info[["Sample Information"]] <- NULL
        debug_info[[".Total IDs"]] <- sprintf("Found: %d", length(sample_ids))
        if (length(sample_ids) > 0) {
            for (i in seq_along(head(sample_ids, 3))) {
                debug_info[[sprintf("..Sample ID %d", i)]] <- sample_ids[i]
            }
        }

        # Print debug information using the new function
        print_debug_info(debug_info)
    }

#    
#    metadata <- load_and_process_experiment_metadata(
#        metadata_path = metadata_path,
#        categories = EXPERIMENT_CONFIG$CATEGORIES,
#        column_order = EXPERIMENT_CONFIG$COLUMN_ORDER,
#        sample_ids = sample_ids,
#        stop_on_missing_columns = TRUE
#    )
#    
#    # Create color scheme
#    color_scheme <- create_color_scheme(
#        config = list(
#            placeholder = GENOME_TRACK_CONFIG$color_placeholder,
#            input = GENOME_TRACK_CONFIG$color_input
#        ),
#        categories = list(
#            antibody = unique(metadata$antibody),
#            rescue_allele = unique(metadata$rescue_allele)
#        )
#    )
#    
#    # Add after processing metadata but before track creation
#    short_sample_ids <- create_minimal_identifiers(
#        metadata$sample_id,
#        verbose = RUNTIME_CONFIG$debug_verbose
#    )
#    
#    # Create mapping between full and short IDs
#    sample_id_mapping <- setNames(
#        short_sample_ids,
#        metadata$sample_id
#    )
#
}
