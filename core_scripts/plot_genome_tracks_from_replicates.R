#!/usr/bin/env Rscript
################################################################################
# Plot bigwig files that are replicates
# Author: Luis | Date: 2025-05-02 | Version: 2.0.0
################################################################################
# PURPOSE: Plot replicates from multiple experiment-ids
#
# USAGE:
#   ./core_scripts/plot_genome_tracks_from_replicates.R --experiment-id=<experiment-id>
#   experiment-ids should be provided as csv
#
# DEPENDENCIES:
#   configuration_experiment_bmc
#   ~/lab_utils/core_scripts/setup_bmc_experiment.R outputs
#   ~/lab_utils/core_scripts/{logging,script_control,file_operations}.R
#   required_packages
#
# OUTPUTS:
# - Svg or pdf files with genome tracks based on replicate information in the configuration_experiment_bmc and experiment-ids user provides.
################################################################################

#---------------------------------------
# Load user functions
#---------------------------------------
function_filenames <- c("logging", "script_control", "file_operations")
for (function_filename in function_filenames) {
    function_filepath <- sprintf("~/lab_utils/core_scripts/functions_for_%s.R", function_filename)
    normalized_path <- normalizePath(function_filepath)
    if (!file.exists(normalized_path)) {
        stop(sprintf("[FATAL] File with functions not found: %s", normalized_path))
    }
    source(normalized_path)
}

# Proceed if packages are installed. Can be disable.
required_packages <- c("rtracklayer", "GenomicRanges", "Gviz")
check_required_packages(required_packages, verbose = TRUE, skip_validation = args$skip_validation)

#-------------------------------------------------------------------------------
# Handle script arguments
#-------------------------------------------------------------------------------
# Parse arguments and validate configurations
description <- "Plot tracks from different experiments in same plot."
args <- parse_common_arguments(description = description)
experiment_id <- args$experiment_id
accept_configuration <- args$accept_configuration
experiment_dir <- args$experiment_dir
is_template <- args$is_template
file_directory <- if (is_template) args$experiment_dir else file.path(args$experiment_dir, "documentation")
file_identifier <- if (is_template) "template" else args$experiment_id
config_path <- file.path(file_directory, paste0(file_identifier, "_configuration_experiment_bmc"))
metadata_path <- file.path(file_directory, paste0(file_identifier, "_sample_grid.csv"))
args_info <- list(
    title = "Script Configuration",
    "script.name" = get_script_name(),
    "script.description" = description
)
print_debug_info(modifyList(args_info, args))

#-------------------------------------------------------------------------------
# Load and Validate Experiment Configuration and Dependencies
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
# Setup genome and feature files
#-------------------------------------------------------------------------------
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

#--------------------
# Load the metadata dataframe with all experiments.
#--------------------
project_id <- c()
project_bigwig_files <- c()
project_metadata <- data.frame()
total_expected_rows <- 0
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
    if (length(unique(project_id)) != 1) {
        warning("All experiments must belong to the same project. Found projects: ",
             paste(unique(project_id), collapse = ", "))
    }
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

    metadata <- load_and_process_experiment_metadata(
        metadata_path = metadata_path,
        categories = EXPERIMENT_CONFIG$CATEGORIES,
        column_order = EXPERIMENT_CONFIG$COLUMN_ORDER,
        sample_ids = sample_ids,
        stop_on_missing_columns = TRUE
    )

    total_expected_rows <- total_expected_rows + nrow(metadata)
    metadata$experiment_id <- EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID
    project_metadata <- rbind(project_metadata, metadata)
    project_bigwig_files <- c(project_bigwig_files, bigwig_files)

} # end for loop
output_dir <- file.path(Sys.getenv("HOME"), "data", unique(project_id), "plots", "genome_tracks", "replicates")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (nrow(project_metadata) == 0) {
    stop("No metadata entries found after merging all experiments")
}

if (any(duplicated(project_metadata$sample_id))) {
    warning("Duplicate sample IDs found in merged project_metadata")
}

if (nrow(project_metadata) != total_expected_rows) {
    stop(sprintf(
        "project_metadata merge error: Expected %d rows but got %d rows",
        total_expected_rows,
        nrow(project_metadata)
    ))
}

if (length(unique(project_metadata$experiment_id)) != length(config_path)) {
    warning("Number of unique experiment IDs doesn't match number of config files processed")
}

if(RUNTIME_CONFIG$debug_verbose) {
    debug_info <- list(
        "title" = "Metadata Summary",

        "Dataset Metrics" = NULL,
        ".Samples" = sprintf("%d samples, %d columns", nrow(project_metadata), ncol(project_metadata)),
        ".Experiment" = unique(project_metadata$experiment_id),

        "Sample Categories" = NULL,
        ".Types" = paste(
            sprintf("%s: %d", 
            names(table(project_metadata$sample_type)),
            table(project_metadata$sample_type)
            ),
            collapse = ", "
        ),
        ".Antibodies" = paste(
            sprintf("%s: %d", 
                   names(table(project_metadata$antibody)), 
                   table(project_metadata$antibody)
            ), 
            collapse = ", "
        ),
        "Experimental Conditions" = NULL,
        ".Rescue Alleles" = paste(
            sprintf("%s: %d", 
                   names(table(project_metadata$rescue_allele)), 
                   table(project_metadata$rescue_allele)
            ), 
            collapse = ", "
        ),
        ".Treatments" = paste(
            sprintf("%s: %d", 
            names(table(project_metadata$auxin_treatment)),
            table(project_metadata$auxin_treatment)
            ),
            collapse = ", "
        ),

        "ID Validation" = NULL,
        ".Sample ID Range" = sprintf("%s - %s", min(project_metadata$sample_id), max(project_metadata$sample_id)),
        ".Name Format" = sprintf("Example: %s  %s", project_metadata$full_name[1], project_metadata$short_name[1])
    )
    print_debug_info(debug_info)

}

# Create color scheme
color_scheme <- create_color_scheme(
    config = list(
        placeholder = GENOME_TRACK_CONFIG$color_placeholder,
        input = GENOME_TRACK_CONFIG$color_input
    ),
    categories = list(
        antibody = unique(project_metadata$antibody),
        rescue_allele = unique(project_metadata$rescue_allele)
    )
)

# Add after processing metadata but before track creation
short_sample_ids <- create_minimal_identifiers(
    project_metadata$sample_id,
    verbose = RUNTIME_CONFIG$debug_verbose
)

# Create mapping between full and short IDs
sample_id_mapping <- setNames(
    short_sample_ids,
    project_metadata$sample_id
)

if (!is.null(args$override)) {
    structured_log_info("Reoverriding for sample processing")
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

# Determine which groups to process
unique_short_names <- unique(project_metadata$short_name)

if (length(unique_short_names) == 0) {
    stop("No short names found in metadata")
}

# Validate process_file_index if single file processing
#if (RUNTIME_CONFIG$process_single_file && 
#    (RUNTIME_CONFIG$process_file_index <= 0 || 
#     RUNTIME_CONFIG$process_file_index > length(unique_short_names))) {
#    stop(sprintf(
#        "Invalid process_file_index: %d. Must be between 1 and %d",
#        RUNTIME_CONFIG$process_file_index,
#        length(unique_short_names)
#    ))
#}

# Determine processing indices
short_name_to_process <- if (RUNTIME_CONFIG$process_single_file) {
    RUNTIME_CONFIG$process_file_index
} else {
    seq_along(unique_short_names)
}

scaling_modes <- c("local", "individual")
# Process each group
for (short_name_idx in short_name_to_process) {
    current_short_name <- unique_short_names[short_name_idx]

    # Use exact matching with the short name
    row_samples_to_visualize <- project_metadata[project_metadata$short_name == current_short_name, ]
    row_samples_to_visualize <- row_samples_to_visualize[order(row_samples_to_visualize$experiment_id), ]

    if (nrow(row_samples_to_visualize) == 0) {
        warning(sprintf("No samples found for short name: %s", current_short_name))
        next
    }

    if (RUNTIME_CONFIG$debug_verbose) {
        message(sprintf(
            "Processing group %d/%d: %s (%d samples)",
            short_name_idx,
            length(unique_short_names),
            current_short_name,
            nrow(row_samples_to_visualize)
        ))
    }

    label_result <- create_track_labels(
        samples = row_samples_to_visualize,
        always_show = GENOME_TRACK_CONFIG$label_always_show,
        never_show = GENOME_TRACK_CONFIG$label_never_show,
        separator = GENOME_TRACK_CONFIG$label_separator,
        verbose = RUNTIME_CONFIG$debug_verbose
    )

    if (!label_result$success) {
        warning("Failed to create track labels for group ", short_name_idx, ": ", label_result$error)
        track_labels <- row_samples_to_visualize$short_name
        if (RUNTIME_CONFIG$debug_verbose) {
            message("  Using fallback labels (short_name)")
        }
    } else {
        track_labels <- label_result$data$labels
        if (RUNTIME_CONFIG$debug_verbose) {
            message("\nTrack Labels Created:")
            invisible(mapply(function(id, label) {
                message(sprintf("  %s -> %s", id, label))
            }, row_samples_to_visualize$sample_id, track_labels))
        }
    }

    if (length(track_labels) != nrow(row_samples_to_visualize)){
        warning("number of track labels does not match number of row samples")
    }

    # Initialize tracks list with chromosome axis
    tracks <- list(
        Gviz::GenomeAxisTrack(
            name = sprintf(GENOME_TRACK_CONFIG$format_genome_axis_track_name, chromosome_to_plot)
        )
    )

    for (i in seq_len(nrow(row_samples_to_visualize))) {
        sample_id <- row_samples_to_visualize$sample_id[i]
        current_antibody <- row_samples_to_visualize$antibody[i]
        track_label <- track_labels[i]
        bigwig_file_path <- project_bigwig_files[grepl(sample_id, project_bigwig_files)][1]
        if (RUNTIME_CONFIG$debug_verbose) {
            debug_info <- list(
                "title" = "Sample Processing and Track Creation",

                "Processing Status" = NULL,
                ".Progress" = sprintf("Sample %d/%d", i, nrow(row_samples_to_visualize)),
                ".Sample ID" = sample_id,
                ".Antibody" = current_antibody,
                ".Track Label" = track_label,

                "Bigwig File Matching" = NULL,
                ".Available Files" = sprintf("Total: %d", length(project_bigwig_files))
            )

            # Add first few available files
            for (i in seq_along(head(project_bigwig_files, 3))) {
                debug_info[[sprintf("..File %d", i)]] <- basename(project_bigwig_files[i])
            }

            # Add matching results
            matches <- grepl(sample_id, project_bigwig_files)
            debug_info[[".Pattern Matches"]] <- sprintf("Found: %d", sum(matches))

            # Add matching files if any
            if (sum(matches) > 0) {
                debug_info[[".Matching Files"]] <- NULL
                for (i in seq_along(head(project_bigwig_files[matches], 3))) {
                    debug_info[[sprintf("..Match %d", i)]] <- basename(project_bigwig_files[matches][i])
                }
            }

            # Track creation details
            debug_info[["Track Configuration"]] <- NULL
            debug_info[[".Bigwig Status"]] <- if(is.na(bigwig_file_path)) "NOT FOUND" else basename(bigwig_file_path)
            debug_info[[".Sample Mapping"]] <- sample_id_mapping[sample_id]
            debug_info[[". Second Verification Track Label"]] <- track_label

            print_debug_info(debug_info)
        }

        if (length(bigwig_file_path) != 1) {
            message(sprintf("Bigwig files found does not equal 1: %s", current_short_name))
            warning(sprintf("Sample id used", sample_id))
            next
        }

        track_name_arguments <- c(
            sample_id_mapping[sample_id],  # Mapped ID
            track_label                # Generated label from create_track_labels
        )

        placeholder_name_arguments <- c(
            sample_id_mapping[sample_id],
            track_label,
            GENOME_TRACK_CONFIG$format_suffix
        )

        track_color <- if (current_antibody == "Input") {
            color_scheme$fixed$input
        } else {
            color_scheme$get_color("antibody", current_antibody)
        }

        track_creation_result <- create_sample_track(
            bigwig_file_path = bigwig_file_path,
            track_format_name = GENOME_TRACK_CONFIG$format_sample_track_name,
            format_args = track_name_arguments,
            track_color = track_color,
            track_type = GENOME_TRACK_CONFIG$track_defaults_sample$type,
            genomic_range = genome_range,
            track_params = GENOME_TRACK_CONFIG$track_defaults_sample,
            verbose = RUNTIME_CONFIG$debug_verbose
        )

        if (track_creation_result$success) {
            if (RUNTIME_CONFIG$debug_verbose) {
                message(sprintf("  Successfully created track for sample: %s", sample_id))
            }
            tracks[[length(tracks) + 1]] <- track_creation_result$data
        } else {
            # Report the error
            message(sprintf("\nTrack creation failed for sample %s:", sample_id))
            message(sprintf("  Error: %s", track_creation_result$error))
            message("  Creating placeholder track instead")
            if (RUNTIME_CONFIG$debug_verbose) {
                # Add context about the failure
                message("  Context:")
                message(sprintf("    Genomic Range: %s", 
                               if(exists("genome_range")) "Present" else "Missing"))
                message(sprintf("    Chromosome: %s", chromosome_roman))
                message(sprintf("    Width: %d", chromosome_width))
            }
            # Create placeholder with error handling
            placeholder_track_creation_result <- create_placeholder_track(
                sampling_rate = GENOME_TRACK_CONFIG$track_sampling_rate,
                chromosome_width = chromosome_width,
                track_color = GENOME_TRACK_CONFIG$color_placeholder,
                type = GENOME_TRACK_CONFIG$track_type,
                chromosome_name = chromosome_roman,
                placeholder_format_name = GENOME_TRACK_CONFIG$format_placeholder_track_name,
                format_args = placeholder_name_arguments,
                track_params = GENOME_TRACK_CONFIG$track_defaults_placeholder,
                verbose = RUNTIME_CONFIG$debug_verbose
            )

            if (!placeholder_track_creation_result$success) {
                stop(sprintf("Failed to create placeholder track: %s", placeholder_track_creation_result$error))
            } else {
                tracks[[length(tracks) + 1]] <- placeholder_track_creation_result$data
            }
        }
    }

    if (exists("reference_grange")) {
        if (!is(reference_grange, "GRanges")) {
            stop("reference_grange must be a GRanges object")
        }

        track_args <- c(
            list(
                range = reference_grange,
                name = "Reference"
            ),
            GENOME_TRACK_CONFIG$track_defaults_sample
        )

        # Optional: Debug info
        if (RUNTIME_CONFIG$debug_verbose) {
            message("Creating reference track with arguments:")
            #print(str(track_args))
        }

        tracks[[length(tracks) + 1]] <- do.call(
            Gviz::DataTrack,
            track_args
        )
    }

    if (exists("features")) {
        tracks[[length(tracks) + 1]] <- Gviz::AnnotationTrack(
            features,
            name = "Features",
            size = 0.5,
            background.title = "lightgray",
            fontcolor.title = "black",
            cex.title = 0.6
        )
    }

    title_replicate_template <- paste(
        "Project_id: %s",
        "Replicate: %s",
        "Chromsome: %s",
        "Time: %s",
        "Normalization:%s",
        sep = "\n"

    )

    plot_title <- sprintf(
        title_replicate_template,
        unique(project_id),
        current_short_name,
        chromosome_to_plot,
        TIME_CONFIG$current_timestamp,
        normalization_method
    )

    plot_config <- create_track_plot_config(
        tracks = tracks,
        chromosome = chromosome_roman,
        to = chromosome_width,
        title = plot_title,
        visualization_params = GENOME_TRACK_CONFIG$plot_defaults,
        verbose = RUNTIME_CONFIG$debug_verbose
    )

    for (mode in scaling_modes) {
        # Create filename for this mode
        filename_format_comparison_templates <- "%s_%s_idx%02d_chr%s_%s.svg"
        plot_filename <- sprintf(
            filename_format_comparison_templates,
            TIME_CONFIG$current_timestamp,
            gsub("_", "", unique(project_id)),
            short_name_idx,
            chromosome_to_plot,
            mode
         )

         plot_file <- file.path(
             output_dir,
             plot_filename
         )

        if (RUNTIME_CONFIG$debug_verbose) {
            message(paste(rep("-", 80), collapse = ""))
            message(sprintf("\nScaling mode: %s", mode))
            message("\nPlot Generation Details:")
            message(sprintf("  Title Components:"))
            message(sprintf("    Experiment: %s", experiment_id))
            message(sprintf("    Chromosome: %s", chromosome_to_plot))
            message(sprintf("    Sample Count: %d", nrow(row_samples_to_visualize)))
            message(sprintf("    Timestamp: %s", TIME_CONFIG$current_timestamp))
            message(sprintf("    Normalization: %s", normalization_method))
            message("\n  Output Configuration:")
            message(sprintf("    Plot Directory: %s", output_dir))
            message(sprintf("    Filename: %s", basename(plot_file)))
            message(sprintf("    Full Path: %s", plot_file))
            # Visual separator for readability in log
            message(paste(rep("-", 80), collapse = ""))
        }

        # Set y-limits based on mode
        if (mode == "local") {
            # Calculate limits for this group's files
            print(row_samples_to_visualize$sample_id)
            #group_files <- project_bigwig_files[grepl(
            #    paste(row_samples_to_visualize$sample_id, collapse = "|"),
            #    project_bigwig_files
            #)]
            pattern <- sprintf("processed_(%s)_sequence_to_S288C_CPM\\.bw$", 
                               paste(row_samples_to_visualize$sample_id, collapse = "|"))
            group_files <- grep(pattern, project_bigwig_files, value = TRUE, perl = TRUE)
            if (RUNTIME_CONFIG$debug_verbose) {
                message("Selected sample IDs: ", paste(row_samples_to_visualize$sample_id, collapse=", "))
                message("Found files: ", paste(basename(group_files), collapse=", "))
                if (length(group_files) != nrow(row_samples_to_visualize)) {
                    warning(sprintf("Number of files (%d) doesn't match number of samples (%d)",
                            length(group_files), nrow(row_samples_to_visualize)))
                }
            }

            local_limits_result <- calculate_track_limits(
                bigwig_files = group_files,
                genome_range = genome_range,
                padding_fraction = .1,
                #padding_fraction = GENOME_TRACK_CONFIG$ylim_padding,
                verbose = RUNTIME_CONFIG$debug_verbose
            )
            plot_config$ylim <- if (local_limits_result$success) 
                local_limits_result$data else GENOME_TRACK_CONFIG$ylim_defaults
        } else {  # individual
            plot_config$ylim <- NULL  # Let tracks scale independently
        }

        if (RUNTIME_CONFIG$output_dry_run) {
            # Display only
            execute_track_plot(
                plot_config = plot_config,
                plot_params = GENOME_TRACK_CONFIG$plot_defaults,
                display_plot = TRUE,
                verbose = RUNTIME_CONFIG$debug_verbose
            )
        } else {
            # Save plot
            execute_track_plot(
                plot_config = plot_config,
                save_path = plot_file,
                save_params = list(
                    width = GENOME_TRACK_CONFIG$display_width,
                    height = GENOME_TRACK_CONFIG$display_height
                ),
                plot_params = GENOME_TRACK_CONFIG$plot_defaults,
                display_plot = FALSE,
                verbose = RUNTIME_CONFIG$debug_verbose
            )
        }
    }
}
