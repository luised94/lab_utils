################################################################################
# Parse FastQC files for ChIP-seq quality control
################################################################################
# PURPOSE:
# USAGE:
# !! ----> REQUIRED UPDATES:
# STRUCTURE:
# VALIDATION:
# DEPENDENCIES:
# COMMON ISSUES:
# OUTPUT:
# AUTHOR: Luis
# DATE: 
# VERSION: 1.0.0
################################################################################
################################################################################
# Configuration and Debug Settings
################################################################################
DEBUG_CONFIG <- list( # !! UPDATE THIS
    single_file_mode = FALSE,           # Test single file in main logic.
    verbose = TRUE,           # Print processing details
    interactive = TRUE,       # Allow interactive processing
    dry_run = FALSE,         # Skip file writes
    files_to_process_idx = 1  # Process specific files in debug mode
)

SCRIPT_SPECIFIC_CONFIG <- list(

)
TIME_CONFIG <- list(
    timestamp_format = "%Y%m%d_%H%M%S",
    date_format = "%Y%m%d"
)

# Generate timestamps
TIMESTAMPS <- list(
    full = format(Sys.time(), TIME_CONFIG$timestamp_format),
    date = format(Sys.Date(), TIME_CONFIG$date_format)
)
################################################################################
# Load and Validate Experiment Configuration
################################################################################
# Bootstrap phase
bootstrap_path <- normalizePath("~/lab_utils/core_scripts/functions_for_file_operations.R", 
                              mustWork = FALSE)
if (!file.exists(bootstrap_path)) {
    stop(sprintf("[FATAL] Bootstrap file not found: %s", bootstrap_path))
}
source(bootstrap_path)
# Define required dependencies
required_modules <- list(
    list(
        path = "~/lab_utils/failsafe_scripts/functions_for_logging.R",
        description = "BMC Configuration",
        required = TRUE
    )
)
################################################################################
# Directory Setup and Validation
################################################################################
experiment_id <- "241007Bel"  # !! UPDATE THIS
base_dir <- file.path(Sys.getenv("HOME"), "data", experiment_id)
qc_dir <- file.path(base_dir, FASTQC_CONFIG$qc_subdir)

stopifnot(
    "Base directory does not exist" = dir.exists(base_dir),
    "Quality control directory does not exist" = dir.exists(qc_dir)
)
################################################################################
# File Discovery
################################################################################
fastqc_files <- list.files(
    qc_dir,
    pattern = FASTQC_CONFIG$fastqc_pattern,
    recursive = TRUE,
    full.names = TRUE
)

stopifnot(
    "No FastQC files found" = length(fastqc_files) > 0
)
################################################################################
# Process Files
################################################################################
files_to_process <- if (DEBUG_CONFIG$single_file_mode) {
    DEBUG_CONFIG$files_to_process_idx
} else {
    seq_along(fastqc_files)
}

stopifnot(
    `Files to process must be within valid range` = all(files_to_process <= length(fastqc_files)),
    `Files to process cannot be negative` = all(files_to_process > 0)
)

for (file in files_to_analyze){
    # specific analysis code.
}
