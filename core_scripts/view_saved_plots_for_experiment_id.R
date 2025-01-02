#!/usr/bin/env Rscript
# Required Packages and Functions
source(file.path(Sys.getenv("HOME"), "lab_utils", "core_scripts", "functions_for_script_control.R"))

# Parse arguments and validate configurations
args <- parse_args(commandArgs(trailingOnly = TRUE))
experiment_id <- args[["experiment-id"]]
source(file.path("~/data", experiment_id, "documentation", 
                paste0(experiment_id, "_bmc_config.R")))
validate_configs(c("RUNTIME_CONFIG", "EXPERIMENT_CONFIG"))
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
