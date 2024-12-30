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
    base_dir = VIEWER_CONFIG$base_dir,
    patterns = VIEWER_CONFIG$patterns,
    experiment = "241007Bel",
    #timestamp = "20231116",  # Optional
    #pattern = "chr10",      # Optional
    verbose = RUNTIME_CONFIG$verbose
)

if (length(files) > 0) {
    if (RUNTIME_CONFIG$verbose) {
        base::message("\nStarting plot display...")
    }
    
    # Usage in main script
    display_plots(
        files = files,
        device_config = VIEWER_CONFIG$device,
        interactive = RUNTIME_CONFIG$interactive,
        display_time = RUNTIME_CONFIG$display_time,
        verbose = RUNTIME_CONFIG$verbose
    )
} else {
    base::message("No files found matching criteria")
}
