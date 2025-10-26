#!/usr/bin/env Rscript
if (!interactive()) {
    message("\nThis script is designed to be run interactively.")
    message("\nTo use this script:")
    message("1. Start R in interactive mode:")
    message("   R")
    message("\n2. Set required variables:")
    message('   experiment_id <- "241010Bel"')
    message("   chromosome_to_plot <- 10")
    message("   samples_per_batch <- 4")
    message("\n3. Source the script:")
    message('   source("plot_genome_tracks_in_batches.R")')
    message("\nAlternatively, use the companion script for batch processing:")
    message('   Rscript plot_genome_tracks_batch_mode.R --experiment-id "241010Bel" ...')
    quit(status = 1)
}
## Required Packages and Functions
source(file.path(Sys.getenv("HOME"), "lab_utils", "core_scripts", "functions_for_script_control.R"))

# At script start
required_variables <- c(
    "experiment_id",
    "chromosome_to_plot",
    "samples_per_batch"
)

# Check and guide
missing_vars <- required_variables[!sapply(required_variables, exists)]
if (length(missing_vars) > 0) {
    message("\nMissing required variables:")
    message(paste(" -", missing_vars, collapse = "\n"))
    message("\nPlease set these variables and re-source the script:")
    message("Example:")
    message('experiment_id <- "241010Bel"')
    message("chromosome_to_plot <- 10")
    message("samples_per_batch <- 4")
    message("\nThen:")
    message('source("plot_genome_tracks_in_batches.R")')
    return(invisible(NULL))
}

# Comprehensive variable validation
stopifnot(
    # experiment_id validation
    "experiment_id must be character" = 
        is.character(experiment_id),
    "experiment_id must be a single value" = 
        length(experiment_id) == 1,
    "experiment_id must match pattern YYMMDD'Bel'" = 
        grepl("^\\d{6}Bel$", experiment_id),
    
    # chromosome_to_plot validation
    "chromosome_to_plot must be numeric" = 
        is.numeric(chromosome_to_plot),
    "chromosome_to_plot must be a single value" = 
        length(chromosome_to_plot) == 1,
    "chromosome_to_plot must be between 1 and 16" = 
        chromosome_to_plot >= 1 && chromosome_to_plot <= 16,
    
    # samples_per_batch validation
    "samples_per_batch must be numeric" = 
        is.numeric(samples_per_batch),
    "samples_per_batch must be a single value" = 
        length(samples_per_batch) == 1,
    "samples_per_batch must be positive" = 
        samples_per_batch > 0,
    "samples_per_batch must be an integer" = 
        samples_per_batch == round(samples_per_batch),
    "samples_per_batch must be reasonable (1-20)" = 
        samples_per_batch >= 1 && samples_per_batch <= 20
)

# If we get here, all validations passed
message("\nVariable validation successful:")
message(sprintf("  Experiment: %s", experiment_id))
message(sprintf("  Chromosome: %d", chromosome_to_plot))
message(sprintf("  Samples per batch: %d", samples_per_batch))
#
## Parse arguments and validate configurations
#args <- parse_args(commandArgs(trailingOnly = TRUE))
#experiment_id <- args[["experiment-id"]]
#source(file.path("~/data", experiment_id, "documentation", 
#                paste0(experiment_id, "_configuration_experiment_bmc")))
#validate_configs(c("RUNTIME_CONFIG", "EXPERIMENT_CONFIG"))
#-----------------------------------------------------------------------------
required_packages <- c("rsvg", "magick")
for (pkg in required_packages) {
    if (!base::requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' is missing", pkg))
    }
}

source("~/lab_utils/core_scripts/functions_for_plotting_utilities.R")
# Main Processing
#-----------------------------------------------------------------------------
# Command line argument processing could be added here

files <- find_plot_files(
    base_dir = VIEWER_CONFIG$path_base,
    patterns = VIEWER_CONFIG$pattern_svg,
    experiment = "241007Bel",
    #timestamp = "20231116",  # Optional
    #pattern = "chr10",      # Optional
    verbose = RUNTIME_CONFIG$debug_verbose
)

if (length(files) > 0) {
    if (RUNTIME_CONFIG$debug_verbose) {
        base::message("\nStarting plot display...")
    }
    
    # Usage in main script
    display_plots(
        files = files,
        device_config = list(
            width = VIEWER_CONFIG$display_width,
            height = VIEWER_CONFIG$display_height,
        ),
        interactive = RUNTIME_CONFIG$debug_interactive,
        display_time = RUNTIME_CONFIG$output_display_time,
        verbose = RUNTIME_CONFIG$debug_verbose
    )
} else {
    base::message("No files found matching criteria")
}
