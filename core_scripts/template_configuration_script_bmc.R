################################################################################
# Ineractive script configuration for BMC experiments
# Author: Luis | Date: 2025-05 | Version: 2.0.0
################################################################################
# PURPOSE:
#   Contains configuration and parameters required for interactive scripts to analyze bmc experiments
#
# USAGE:
# 1. Update experiment_id (format: YYMMDDBel, e.g., "241122Bel")
# 2. Update other key parameters.
# 3. Source script.
#
# DEPENDENCIES: 
#   ~/lab_utils/core_scripts/override_configuration.R
#
# OUTPUTS:
#   Depends on the script.
################################################################################
message("Sourcing script configuration file...")
#---------------------
# CONFIGURATION_PARAMETERS
#---------------------
# --- Variable ---
# Likely updated every time a new experiment is created.
# Write single experiment id or single character that is comma-separated
EXPERIMENT_IDS <- "250324Bel"
EXPERIMENT_ID <- "Exp_20250515_1"

target_comparison_columns <- c("rescue_allele", "suppressor_allele")
columns_to_exclude_from_replicate_determination <- c(
  "sample_type", "sample_ids",
  "bigwig_file_paths", "full_name",
  "short_name", "experiment_id",
  "repeats"
)

ROOT_DIRECTORY <- system("git rev-parse --show-toplevel", intern = TRUE)
#system("git rev-parse --is-inside-work-tree", intern = TRUE)
#if command -v git &>/dev/null && git rev-parse --is-inside-work-tree &>/dev/null; then
# --- Constant ---
#  git rev-parse --show-toplevel
#  return 0
#fi
# TODO: Need to do global but not configured yet
# TODO: Remove and rework these settings.
ACCEPT_CONFIGURATION <- TRUE
BAM_PROCESSING <- "blFiltered"
BIGWIG_NORM_METHOD <- "RAW"
CHROMOSOMES_TO_PLOT <- c(7, 10, 14)
EXPECTED_FORMAT_EXPERIMENT_ID <- "^\\d{6}Bel$"
FASTQ_PATTERN <- "consolidated_.*_sequence\\.fastq$"
GENOME_TRACK_Y_AXIS_SCALING <- c("individual")
OUTPUT_FORMAT <- "svg"
OVERRIDE_CONFIGURATION_PATH <- "~/lab_utils/core_scripts/override_configuration.R"
PADDING_FRACTION <- 0.1
SAMPLE_ID_CAPTURE_PATTERN <- "consolidated_([0-9]{1,6})_sequence\\.fastq$"
SKIP_PACKAGE_CHECKS <- TRUE
VALID_GENOME_TRACK_SCALING_MODES <- c("local", "individual")
VALID_OUTPUT_FORMATS <- c("svg", "pdf", "png")
VARIABLES_TO_REMOVE <- c("IS_COMMA_SEPARATED", "missing_dirs")

#TIME_CONFIG <- list(
#  # Format specifications
#  timestamp_format = "%Y%m%d_%H%M%S",  # YYYYMMDD_HHMMSS
#  date_format = "%Y%m%d",        # YYYYMMDD
#
#  # Current values
#  current_timestamp = format(Sys.time(), "%Y%m%d_%H%M%S"),
#  current_date = format(Sys.Date(), "%Y%m%d")
#)

#EXPERIMENT_DIR <-
#LABEL_MAPPINGS <- list()
#REQUIRED_DIRECTORIES <-
##LOG_TO_FILE <-
#OUTPUT_DIR <-
#OVERRIDE <- NULL
#SCRIPT_TO_RUN <- ""
#---------------------
# Validation layer
#---------------------
stopifnot(
  "EXPERIMENT_IDS is required" =
    !is.null(EXPERIMENT_IDS),
  "EXPERIMENT_IDS should be character vector of length 1." =
    length(EXPERIMENT_IDS) == 1,
  "OUTPUT_FORMAT must be svg, pdf or png." =
    OUTPUT_FORMAT %in% VALID_OUTPUT_FORMATS
)

#---------------------
# EXPERIMENT CONFIGURATION SETUP
#---------------------
OUTPUT_EXTENSION <- paste0(".", OUTPUT_FORMAT)
BIGWIG_PATTERN <- sprintf(
  fmt = "processed_.*_sequence_to_S288C_%s_%s\\.bw$",
  BAM_PROCESSING, BIGWIG_NORM_METHOD
)
reproducible_subset_quote_list <- "~/lab_utils/core_scripts/metadata_subset.R"
FILE_GENOME_DIRECTORY <- file.path(Sys.getenv("HOME"), "data", "REFGENS")
FILE_GENOME_PATTERN <- "S288C_refgenome.fna"
FILE_FEATURE_DIRECTORY <- file.path(Sys.getenv("HOME"), "data", "feature_files")
FILE_FEATURE_PATTERN <- "eaton_peaks"
row_filtering_expression <- quote(!(rescue_allele == "4R" & suppressor_allele == "NONE"))

IS_COMMA_SEPARATED <- grepl(",", EXPERIMENT_IDS)
# Setup experiment directories ----------------
if (!IS_COMMA_SEPARATED){
  EXPERIMENT_DIR <- normalizePath(
    file.path(Sys.getenv("HOME"), "data", EXPERIMENT_IDS),
    mustWork = FALSE
  )
}

if (IS_COMMA_SEPARATED) {
  # Split and clean the IDs
  split_experiment_ids <- stri_split_fixed(EXPERIMENT_IDS, ",")[[1]]
  clean_experiment_ids <- trimws(split_experiment_ids)  # Remove any whitespace
  EXPERIMENT_IDS <-  clean_experiment_ids[clean_experiment_ids != ""]  # Remove empty elements
  # Check for duplicates
  if (length(unique(EXPERIMENT_IDS)) != length(EXPERIMENT_IDS)) {
    stop("Duplicate experiment IDs detected")
  }
  # Validate format of each ID
  invalid_ids <- EXPERIMENT_IDS[!grepl(
    EXPECTED_FORMAT_EXPERIMENT_ID,
    EXPERIMENT_IDS,
    perl = TRUE
  )]
  if (length(invalid_ids) > 0) {
    stop(sprintf(
      "Invalid experiment-id format(s):\n%s\nExpected format: YYMMDD'Bel'",
      paste(invalid_ids, collapse = ", ")
    ))
  }
  EXPERIMENT_IDS <- EXPERIMENT_IDS
  EXPERIMENT_DIR <- sapply(EXPERIMENT_IDS, function(experiment_id) {
    normalizePath(
      file.path(Sys.getenv("HOME"), "data", experiment_id),
      mustWork = FALSE
    )
  })
  VARIABLES_TO_REMOVE <- c(VARIABLES_TO_REMOVE,
    "split_experiment_ids", "clean_experiment_ids",
    "invalid_ids")
} # end if comma processing statement

if (!all(grepl(EXPECTED_FORMAT_EXPERIMENT_ID, EXPERIMENT_IDS, perl = TRUE))){
  stop(sprintf(
    fmt = "Invalid experiment-id format.\nExpected: YYMMDD'Bel' or '<exp_id1>,<exp_id2>,...'\nReceived: %s",
    EXPERIMENT_IDS
  ))
}

# Identify missing experiment directories
missing_dirs <- EXPERIMENT_DIR[!dir.exists(EXPERIMENT_DIR)]
if ( length(missing_dirs) > 0 ) {
  warning("The following directories are missing:\n",
        paste(missing_dirs, collapse = "\n")
  )
}

message("All experiment directories exist...")
# end experiment directory setup -------------------
#---------------------
# Clean up
#---------------------
rm(list = VARIABLES_TO_REMOVE)
message("Additional variables removed...")
# end clean up section -------------------

message("Configuration complete...")
