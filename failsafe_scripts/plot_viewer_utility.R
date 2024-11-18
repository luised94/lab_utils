#!/usr/bin/env Rscript

# Plot Viewer Configuration
#-----------------------------------------------------------------------------
DEBUG_CONFIG <- list(
    enabled = FALSE,          # Enable debug output
    verbose = TRUE,           # Print processing details
    display_time = 2,         # Seconds between auto-advance
    interactive = TRUE        # Enable interactive viewing
)

VIEWER_CONFIG <- list(
    base_dir = file.path(Sys.getenv("HOME"), "data"),
    patterns = list(
        svg = "\\.svg$",
        timestamp = "^[0-9]{8}_[0-9]{6}",  # YYYYMMDD_HHMMSS
        experiment = "^[0-9]{6}Bel"
    ),
    device = list(
        width = 10,
        height = 8
    )
)

# Required Packages and Functions
#-----------------------------------------------------------------------------
required_packages <- c("rsvg", "magick")
for (pkg in required_packages) {
    if (!base::requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Package '%s' is missing", pkg))
    }
}

source("~/lab_utils/failsafe_scripts/functions_for_plotting_utilities.R")
# Main Processing
#-----------------------------------------------------------------------------
# Command line argument processing could be added here

# Example usage
files <- find_plot_files(
    base_dir = VIEWER_CONFIG$base_dir,
    experiment = "241010Bel",
    #timestamp = "20231116",  # Optional
    #pattern = "chr10",      # Optional
    verbose = DEBUG_CONFIG$verbose
)

if (length(files) > 0) {
    if (DEBUG_CONFIG$verbose) {
        base::message("\nStarting plot display...")
    }
    
    display_plots(files, VIEWER_CONFIG)
} else {
    base::message("No files found matching criteria")
}
