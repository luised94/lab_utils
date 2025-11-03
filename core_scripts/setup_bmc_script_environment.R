################################################################################
# Centralize entry script for other ngs analysis scripts
# Author: Luis | Date: 2025-11-03 | Version: 1.0.0
################################################################################
# PURPOSE:
#   Sets up execution environment for ngs analysis scripts
#
# USAGE:
# 1. Update configuration files
# 2. Source in other scripts using standard R.
#
# DEPENDENCIES: 
#   ~/lab_utils/core_scripts/override_configuration.R
#
# OUTPUTS:
#   Depends on the script.
################################################################################

if (
  !exists("ROOT_DIRECTORY") ||
  is.null(ROOT_DIRECTORY) ||
  ROOT_DIRECTORY == ""
) {

  ROOT_DIRECTORY <- system("git rev-parse --show-toplevel", intern = TRUE)
    # Check if the command failed (returns error status or empty result)
    if (length(ROOT_DIRECTORY) == 0 || ROOT_DIRECTORY == "") {
      stop("Could not determine git root directory. Not in a git repository? Current working directory: ", getwd())
    }
}


#---------------------
# Validation layer
#---------------------
stopifnot(
  "EXPERIMENT_IDS is required" =
    !is.null(EXPERIMENT_IDS),
  "EXPERIMENT_IDS should be character vector of length 1." =
    length(EXPERIMENT_IDS) == 1,
  "OUTPUT_FORMAT must be svg, pdf or png." =
    OUTPUT_FORMAT %in% OUTPUT_FORMATS_VALID,
  "EXPERIMENT_ID must be in EXPERIMENT_IDS." =
    EXPERIMENT_ID %in% EXPERIMENT_IDS
)

#---------------------
# EXPERIMENT CONFIGURATION SETUP
#---------------------
OUTPUT_EXTENSION <- paste0(".", OUTPUT_FORMAT)
#BIGWIG_PATTERN <- sprintf(
#  fmt = "processed_.*_sequence_to_S288C_%s_%s\\.bw$",
#  BAM_PROCESSING, BIGWIG_NORM_METHOD
#)
BIGWIG_PATTERN <- sprintf(
  "D[0-9]{2}-[0-9]{1,6}_%s\\.bw",
  BIGWIG_NORM_METHOD
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
    file.path(Sys.getenv("HOME"), "data", EXPERIMENT_ID),
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
