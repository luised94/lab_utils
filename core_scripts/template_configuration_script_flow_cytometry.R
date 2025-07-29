################################################################################
# Interactive configuration template for flow cytometry experiments
# Author: Luis | Date: Sometime 2025 | Version: 1.0.0
################################################################################
# PURPOSE:
#   Control settings for flow cytometry analysis scripts.
#
# USAGE:
# 1. Update experiment_id (format: YYMMDDBel, e.g., "241122Bel")
# 2. Source analysis scripts which will reference interactive_flow-cytometry_script_configuration.R
#
# DEPENDENCIES: NONE
#
# OUTPUTS:
# - Output depends on script
################################################################################
message("Sourcing script configuration file...")
#---------------------
# CONFIGURATION_PARAMETERS
#---------------------
# --- Variable ---
# Likely updated every time a new experiment is created.
DIRECTORY_ID <- "250515_G1_arrest_release_arrest"
# Write single experiment id or single character that is comma-separated
EXPERIMENT_ID <- "Exp_20250515_1"

# --- Constant ---
ACCEPT_CONFIGURATION <- TRUE
DROPBOX_PATH <- Sys.getenv("DROPBOX_PATH")
EXPECTED_FORMAT_EXPERIMENT_ID <- "Exp_\\d{8}_\\d{1,6}"
FLOW_CYTOMETRY_BRIDGE_PATH <- "Lab/Experiments/flow_cytometry"
OUTPUT_FORMAT <- "pdf"
OVERRIDE_CONFIGURATION_PATH <- "~/lab_utils/core_scripts/override_configuration.R"
OVERRIDE_THE_CONFIG <- NULL
SKIP_PACKAGE_CHECKS <- TRUE
VALID_OUTPUT_FORMATS <- c("svg", "pdf", "png")
VARIABLES_TO_REMOVE <- c("IS_COMMA_SEPARATED")
#timestamp_format = "%Y%m%d_%H%M%S",    # YYYYMMDD_HHMMSS
#date_format = "%Y%m%d",                # YYYYMMDD
## Current values
#current_timestamp = format(Sys.time(), "%Y%m%d_%H%M%S"),
#current_date = format(Sys.Date(), "%Y%m%d")
#EXPERIMENT_DIR <- 
#LABEL_MAPPINGS <- list()
#REQUIRED_DIRECTORIES <- 
##LOG_TO_FILE <- 
#OUTPUT_DIR <- 
#SCRIPT_TO_RUN <- ""

#---------------------
# Validation layer
#---------------------
stopifnot(
  "EXPERIMENT_ID is required" =
    !is.null(EXPERIMENT_ID),
  "EXPERIMENT_ID should be character vector of length 1." =
    length(EXPERIMENT_ID) == 1,
  "OUTPUT_FORMAT must be svg, pdf or png." =
    OUTPUT_FORMAT %in% VALID_OUTPUT_FORMATS,
  "FLOW_CYTOMETRY_DIR does not exist." =
    dir.exists(FLOW_CYTOMETRY_DIR)
)

if(DROPBOX_PATH == "") {
  message("Environmental variable DROPBOX_PATH not available.")
  message("Either set with my config directory or manually in the parse_flow_cytometry_arguments.")
  stop("!!!! DROPBOX_PATH required for proper directory setting.")
}

#---------------------
# EXPERIMENT CONFIGURATION SETUP
#---------------------
FLOW_CYTOMETRY_DIR <- file.path(DROPBOX_PATH, FLOW_CYTOMETRY_BRIDGE_PATH)
OUTPUT_EXTENSION <- paste0(".", OUTPUT_FORMAT)
IS_COMMA_SEPARATED <- grepl(",", EXPERIMENT_ID)
IS_NOT_VALID_EXPERIMENT_ID <- !all(grepl(EXPECTED_FORMAT_EXPERIMENT_ID, EXPERIMENT_ID, perl = TRUE))

if (!IS_COMMA_SEPARATED){
  SERIES_DIRECTORY <- normalizePath(file.path(FLOW_CYTOMETRY_DIR, DIRECTORY_ID))
}

if (IS_NOT_VALID_EXPERIMENT_ID){
  stop(sprintf(
    fmt = "Invalid experiment-id format.\nExpected: Exp_[0-9]{8}_[0-9] or '<exp_id1>,<exp_id2>,...'\nReceived: %s",
    EXPERIMENT_ID
  ))
}

stopifnot(
  "Series directory does not exists" =
    dir.exists(SERIES_DIRECTORY)
)

message("All experiment directories exist...")

#---------------------
# Clean up
#---------------------
rm(list = VARIABLES_TO_REMOVE)
message("Additional variables removed...")

# Done
message("Configuration complete...")
