#!/usr/bin/env Rscript
# Required Packages and Functions
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
    verbose = DEBUG_CONFIG$verbose
)

if (length(files) > 0) {
    if (DEBUG_CONFIG$verbose) {
        base::message("\nStarting plot display...")
    }
    
    # Usage in main script
    display_plots(
        files = files,
        device_config = VIEWER_CONFIG$device,
        interactive = DEBUG_CONFIG$interactive,
        display_time = DEBUG_CONFIG$display_time,
        verbose = DEBUG_CONFIG$verbose
    )
} else {
    base::message("No files found matching criteria")
}
