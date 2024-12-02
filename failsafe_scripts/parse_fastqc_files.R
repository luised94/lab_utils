################################################################################
# Configuration and Debug Settings
################################################################################
# !! Review debug configuration
DEBUG_CONFIG <- list(
    enabled = FALSE,
    verbose = TRUE,
    interactive = TRUE,
    dry_run = FALSE
)

FASTQC_CONFIG <- list(
    VERSION = "1.0.0",
    FASTQC_DATA_PATTERN = "fastqc_data",
    OUTPUT_SUFFIX = ".tab",
    QC_SUBDIR = "quality_control",
)

#CONFIG$COLUMN_NAMES <- list(
#    SUMMARY = c("Stat", "Value")
#)

# Time formatting configuration
TIME_CONFIG <- list(
    timestamp_format = "%Y%m%d_%H%M%S",  # YYYYMMDD_HHMMSS
    date_format = "%Y%m%d"               # YYYYMMDD
)

# Generate timestamps once at script start
TIMESTAMPS <- list(
    full = format(Sys.time(), TIME_CONFIG$timestamp_format),
    date = format(Sys.Date(), TIME_CONFIG$date_format)
)


experiment_id <- "241122Bel"
base_experiment_dir <- file.path(Sys.getenv("HOME"), "data", experiment_id)
if (is.null(dir_path) || !dir.exists(base_experiment_dir)) {
    stop(sprintf("Invalid directory path: %s", base_experiment_dir))
}
base::sprintf("%s_%s_%s%s",
    get_timestamp(),
    prefix,
    module_name,
    FASTQC_CONFIG$OUTPUT_SUFFIX)


