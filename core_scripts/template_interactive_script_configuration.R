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
EXPECTED_FORMAT_EXPERIMENT_ID <- "^\\d{6}Bel$"

######################
# CONFIGURATION_PARAMETERS
######################
# Write single experiment id or single character that is comma-separated
EXPERIMENT_IDS <- "250324Bel"
CHROMOSOMES_TO_PLOT <- c(7, 10, 14)
VARIABLES_TO_REMOVE <- c("IS_COMMA_SEPARATED", "missing_dirs")
ACCEPT_CONFIGURATION <- TRUE
SKIP_PACKAGE_CHECKS <- TRUE
OUTPUT_FORMAT <- "svg"
OUTPUT_EXTENSION <- paste0(".", OUTPUT_FORMAT)
VALID_OUTPUT_FORMATS <- c("svg", "pdf", "png")
BAM_PROCESSING <- "blFiltered"
BIGWIG_NORM_METHOD <- "RAW"
BIGWIG_PATTERN <- sprintf(
  fmt = "processed_.*_sequence_to_S288C_%s_%s\\.bw$",
  BAM_PROCESSING, BIGWIG_NORM_METHOD
)
FASTQ_PATTERN <- "consolidated_.*_sequence\\.fastq$"
SAMPLE_ID_CAPTURE_PATTERN <- "consolidated_([0-9]{1,6})_sequence\\.fastq$"
GENOME_TRACK_Y_AXIS_SCALING <- c("individual")
# TODO: Need to do global but not configured yet
VALID_GENOME_TRACK_SCALING_MODES <- c("local", "individual")
PADDING_FRACTION <- 0.1
# TODO: Remove and rework these settings.
target_comparison_columns <- c("rescue_allele", "suppressor_allele")
columns_to_exclude_from_replicate_determination <- c(
  "sample_type", "sample_ids",
  "bigwig_file_paths", "full_name",
  "short_name", "experiment_id",
  "repeats"
)
reproducible_subset_quote_list <- "~/lab_utils/core_scripts/metadata_subset.R"
# TODO: Move to configuration. //
FILE_GENOME_DIRECTORY <- file.path(Sys.getenv("HOME"), "data", "REFGENS")
FILE_GENOME_PATTERN <- "S288C_refgenome.fna"
FILE_FEATURE_DIRECTORY <- file.path(Sys.getenv("HOME"), "data", "feature_files")
FILE_FEATURE_PATTERN <- "eaton_peaks"
row_filtering_expression <- quote(!(rescue_allele == "4R" & suppressor_allele == "NONE"))
#EXPERIMENT_DIR <-
#LABEL_MAPPINGS <- list()
#REQUIRED_DIRECTORIES <-
##LOG_TO_FILE <-
#OUTPUT_DIR <-
#OVERRIDE <- NULL
#SCRIPT_TO_RUN <- ""
######################
# Validation layer
######################
stopifnot(
  "EXPERIMENT_IDS is required" =
    !is.null(EXPERIMENT_IDS),
  "EXPERIMENT_IDS should be character vector of length 1." =
    length(EXPERIMENT_IDS) == 1,
  "OUTPUT_FORMAT must be svg, pdf or png." =
    OUTPUT_FORMAT %in% VALID_OUTPUT_FORMATS
)

######################
# EXPERIMENT CONFIGURATION SETUP
######################
IS_COMMA_SEPARATED <- grepl(",", EXPERIMENT_IDS)
# Setup experiment directories ----------------
if (!IS_COMMA_SEPARATED){
  EXPERIMENT_DIR <- normalizePath(file.path(Sys.getenv("HOME"), "data", EXPERIMENT_IDS))
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
    normalizePath(file.path(Sys.getenv("HOME"), "data", experiment_id))
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
  stop("The following directories are missing:\n",
        paste(missing_dirs, collapse = "\n")
  )
}

message("All experiment directories exist...")
# end experiment directory setup -------------------
######################
# Clean up
######################
rm(list = VARIABLES_TO_REMOVE)
message("Additional variables removed...")
# end clean up section -------------------

message("Configuration complete...")
