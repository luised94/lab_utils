# Configuration file for running scripts via source instead of command line
# Substitute arguments provided through cli essentially
# Does not provided arguments supplied by the specific experiments configuration file.
# Scripts that use this file:
#   ./plot_genome_tracks_for_multiple_chromosomes_experiments_and_categories.R
# USAGE: Reference file by using the source function.
message("Sourcing script configuration file...")
######################
# CORE_RESOURCES
######################
OVERRIDE_CONFIGURATION_PATH <- "~/lab_utils/core_scripts/override_configuration.R"

######################
# CORE_PARAMETERS
######################
EXPECTED_FORMAT_EXPERIMENT_ID <- "Exp_\\d{8}_\\d{1,6}"

######################
# CONFIGURATION_PARAMETERS
######################
# Write single experiment id or single character that is comma-separated
EXPERIMENT_ID <- "Exp_20250515_1"
DIRECTORY_ID <- "250303_G1_arrest_degrade_and_release"
ACCEPT_CONFIGURATION <- TRUE
SKIP_PACKAGE_CHECKS <- TRUE
OUTPUT_FORMAT <- "svg"
DROPBOX_PATH <- Sys.getenv("DROPBOX_PATH")
FLOW_CYTOMETRY_BRIDGE_PATH <- "Lab/Experiments/flow_cytometry"
FLOW_CYTOMETRY_DIR <- file.path(DROPBOX_PATH, FLOW_CYTOMETRY_BRIDGE_PATH)
VARIABLES_TO_REMOVE <- c("IS_COMMA_SEPARATED", "missing_dirs")
OVERRIDE_THE_CONFIG <- NULL
#EXPERIMENT_DIR <- 
#LABEL_MAPPINGS <- list()
#REQUIRED_DIRECTORIES <- 
##LOG_TO_FILE <- 
#OUTPUT_DIR <- 
#SCRIPT_TO_RUN <- ""
######################
# Validation layer
######################
stopifnot(
  "EXPERIMENT_ID is required" =
    !is.null(EXPERIMENT_ID),
  "EXPERIMENT_ID should be character vector of length 1." =
    length(EXPERIMENT_ID) == 1,
  "OUTPUT_FORMAT must be svg, pdf or png." =
    OUTPUT_FORMAT %in% c("svg", "pdf", "png"),
  "FLOW_CYTOMETRY_DIR does not exist." =
    dir.exists(FLOW_CYTOMETRY_DIR)
)
if(DROPBOX_PATH == "") {
    message("Environmental variable DROPBOX_PATH not available.")
    message("Either set with my config directory or manually in the parse_flow_cytometry_arguments.")
    stop("!!!! DROPBOX_PATH required for proper directory setting.")
}

######################
# EXPERIMENT CONFIGURATION SETUP
######################
IS_COMMA_SEPARATED <- grepl(",", EXPERIMENT_ID)
# Setup experiment directories ----------------
if (!IS_COMMA_SEPARATED){
  SERIES_DIRECTORY <- normalizePath(file.path(FLOW_CYTOMETRY_DIR, DIRECTORY_ID))
}

#if (IS_COMMA_SEPARATED) {
#  # Split and clean the IDs
#  split_experiment_ids <- stri_split_fixed(EXPERIMENT_ID, ",")[[1]]
#  clean_experiment_ids <- trimws(split_experiment_ids)  # Remove any whitespace
#  EXPERIMENT_ID <-  clean_experiment_ids[clean_experiment_ids != ""]  # Remove empty elements
#  # Check for duplicates
#  if (length(unique(EXPERIMENT_ID)) != length(EXPERIMENT_ID)) {
#    stop("Duplicate experiment IDs detected")
#  }
#  # Validate format of each ID
#  invalid_ids <- EXPERIMENT_ID[!grepl(
#    EXPECTED_FORMAT_EXPERIMENT_ID,
#    EXPERIMENT_ID,
#    perl = TRUE)
#    ]
#  if (length(invalid_ids) > 0) {
#    stop(sprintf(
#      "Invalid experiment-id format(s):\n%s\nExpected format: Exp_[0-9]{8}_[0-9]",
#      paste(invalid_ids, collapse = ", ")
#      ))
#  }
#  EXPERIMENT_ID <- EXPERIMENT_ID
#  EXPERIMENT_DIR <- sapply(EXPERIMENT_ID, function(experiment_id) {
#    normalizePath(file.path(FLOW_CYTOMETRY_DIR, experiment_id))
#  })
#  VARIABLES_TO_REMOVE <- c(VARIABLES_TO_REMOVE,
#    "split_experiment_ids", "clean_experiment_ids",
#    "invalid_ids")
#}

if (!all(grepl(EXPECTED_FORMAT_EXPERIMENT_ID, EXPERIMENT_ID, perl = TRUE))){
  stop(sprintf(
    fmt = "Invalid experiment-id format.\nExpected: Exp_[0-9]{8}_[0-9] or '<exp_id1>,<exp_id2>,...'\nReceived: %s",
    EXPERIMENT_ID
  ))
}


stopifnot(
  "Series directory does not exists" =
    dir.exists(SERIES_DIRECTORY)
)
# Identify missing experiment directories
#missing_dirs <- EXPERIMENT_DIR[!dir.exists(EXPERIMENT_DIR)]
#if (length(missing_dirs) > 0) {
#  # Build a common message for missing directories
#  missing_msg <- sprintf(
#    fmt = "The following directories are missing:\n%s",
#    paste(missing_dirs, collapse = "\n")
#  )
#  # Stop script if directories do not exist.
#  stop(sprintf(
#    fmt = paste("Error: Experiment directories are required for script '%s' to run.\n",
#          "%s\n"),
#    SCRIPT_TO_RUN,
#    missing_msg,
#    ), call. = FALSE
#  )
#  VARIABLES_TO_REMOVE <- c(VARIABLES_TO_REMOVE,
#    "missing_msg")
#}

message("All experiment directories exist...")
# end experiment directory setup -------------------
######################
# Clean up
######################
rm(list = VARIABLES_TO_REMOVE)
message("Additional variables removed...")
# end clean up section -------------------

message("Configuration complete...")
