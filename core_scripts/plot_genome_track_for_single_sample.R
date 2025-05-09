#!/usr/bin/env Rscript
###############################################################################
# Plot bigwig files generated from bam from same sample but differentially processed
################################################################################
# PURPOSE: Plot bigwig files to see if blacklist filtering helps samples with many reads
#          at the ends or in blacklisted regions
# USAGE: source("reference_code/pipeline_completion/plot_bigwig_files.R")
# DEPENDENCIES: GenomicRanges, rtracklayer
# OUTPUT: svg plots with comparisons for cpm/rpkm/raw and for shifted/raw/deduped for cpm counts
# AUTHOR: LEMR
# DATE: 2025-02-25
# DATE_V1_COMPLETE: 2025-04-08
# DATE_V2_COMPLETE: 2025-05-02
################################################################################
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

# Proceed if packages are installed. Can be disable.
required_packages <- c("rtracklayer", "GenomicRanges", "Gviz")
check_required_packages(required_packages, verbose = TRUE, skip_validation = FALSE)

config_path <- "~/lab_utils/core_scripts/template_bmc_config.R"
override_path <- "~/lab_utils/core_scripts/override_configuration.R"
safe_source(config_path)
safe_source(override_path)

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

#todo: Need to fix config requirements for ending with _CONFIG or modify the variable to end with _CONFIG
required_configs <- c("EXPERIMENT_CONFIG", "GENOME_TRACK_CONFIG", "RUNTIME_CONFIG")
validate_configs(required_configs)
invisible(lapply(required_configs, function(config) {
    print_config_settings(get(config), title = config)
}))


# Handle configuration override (independent)
override_result <- apply_runtime_override(
    config = RUNTIME_CONFIG,
    preset_name = "full_inspect_pipeline",
    preset_list = OVERRIDE_PRESETS
)
RUNTIME_CONFIG <- override_result$modified
print_debug_info(modifyList(
    list(
        title = "Final Configuration",
        "override.mode" = override_result$mode
    ),
    RUNTIME_CONFIG
))

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

reference_bigwig <- "~/data/100303Bel/coverage/processed_034475_sequence_to_S288C_RPKM.bw"
stopifnot(
    "Reference bigwig does not exist." = file.exists(reference_bigwig)
)
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


# Initialize tracks list with chromosome axis
tracks <- list(
    Gviz::GenomeAxisTrack(
        name = sprintf(GENOME_TRACK_CONFIG$format_genome_axis_track_name, chromosome_to_plot)
    )
)

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

title_plot_format <- "%s_%s_reference_chr%s_%s.svg"
plot_title <- sprintf(
    title_plot_format,
    TIME_CONFIG$current_timestamp,
    "100303Bel",
    chromosome_to_plot,
    "CPM"
)

plot_config <- create_track_plot_config(
    tracks = tracks,
    chromosome = chromosome_roman,
    to = chromosome_width,
    title = plot_title,
    visualization_params = GENOME_TRACK_CONFIG$plot_defaults,
    verbose = RUNTIME_CONFIG$debug_verbose
)

mode <- "indiviadual"
filename_format_reference_template <- "%s_%s_reference_chr%s_%s.svg"
plot_filename <- sprintf(
    filename_format_reference_template,
    TIME_CONFIG$current_timestamp,
    "100303Bel",
    chromosome_to_plot,
    mode
 )

 plot_file <- file.path(
     "~/",
     plot_filename
 )

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
