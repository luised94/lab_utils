################################################################################
# Ineractive script configuration for BMC experiments
# Author: Luis | Date: 2025-10-27 | Version: 3.0.0
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
#===============================================================================
# GLOBAL CONFIGURATION
#===============================================================================
# Write single experiment id or single character that is comma-separated
EXPERIMENT_IDS <- "250714Bel"
EXPERIMENT_ID <- "250714Bel"

# Execution control
VERBOSE <- TRUE              # Print progress messages to console
DRY_RUN <- FALSE             # Preview operations without executing

# Time and date formats
TIMESTAMP_FORMAT <- "%Y%m%d_%H%M%S"  # YYYYMMDD_HHMMSS
DATE_FORMAT <- "%Y%m%d"              # YYYYMMDD

# TODO: Need to do global but not configured yet
# TODO: Remove and rework these settings.
ACCEPT_CONFIGURATION <- TRUE
BAM_PROCESSING <- "blFiltered"
BIGWIG_NORM_METHOD <- "CPM"
CHROMOSOMES_TO_PLOT <- c(7, 10, 14)
EXPECTED_FORMAT_EXPERIMENT_ID <- "^\\d{6}Bel$"
#FASTQ_PATTERN <- "consolidated_.*_sequence\\.fastq$"
FASTQ_PATTERN <- "fastpfiltered_.*_(R1|R2|NA)\\.fastq$"
#GENOME_TRACK_SCALING_MODES_VALID <- c("local", "individual")
GENOME_TRACK_Y_AXIS_SCALING <- c("individual")
OUTPUT_FORMAT <- "svg"
OUTPUT_FORMATS_VALID <- c("svg", "pdf", "png")
OVERRIDE_CONFIGURATION_PATH <- "~/lab_utils/core_scripts/override_configuration.R"
PADDING_FRACTION <- 0.1
#SAMPLE_ID_CAPTURE_PATTERN <- "consolidated_([0-9]{1,6})_sequence\\.fastq$"
SAMPLE_ID_CAPTURE_PATTERN <- "fastpfiltered_(D[0-9]{2}-[0-9]{1,6})_(R1|R2|NA)\\.fastq$"
SKIP_PACKAGE_CHECKS <- TRUE
VARIABLES_TO_REMOVE <- c("IS_COMMA_SEPARATED", "missing_dirs")
# Define directory structure
DATA_DIRECTORIES <- c(
  "peak",
  "fastq/raw",
  "fastq/processed",
  "quality_control",
  "alignment",
  "coverage",
  "plots/genome_tracks/overview",
  "plots/genome_tracks/experimental_comparisons",
  "documentation/dna_qc_traces",
  "documentation/config"
)

GENOME_TRACK_CONFIG <- list(
  use_custom_visualization = FALSE,  # Control flag

  # Display dimensions
  display_width = 10,
  display_height = 8,

  # Track Creation
  track_points_default = 1000,
  #track_show_title = TRUE,

  # Track defaults
  track_ylim = c(0, 1000),  # Default y-limits, adjust as needed
  track_sampling_rate = 100,  # Points per base pair for empty tracks

  # Track colors
  color_placeholder = "#cccccc",
  color_input = "#808080",

  # Track naming
  format_sample_track_name = "%s: %s",
  format_control_track_name = "%s: %s - %s",
  format_placeholder_track_name = "%s: %s - %s",
  format_suffix = "(No data)",
  format_genome_axis_track_name = "Chr %s Axis",

  # Labels
  label_always_show = "antibody",
  label_never_show = c("sample_id", "full_name", "short_name", "X__cf_genotype"),
  label_separator = "-",

  # File handling
  file_pattern = "consolidated_.*_sequence\\.fastq$",
  file_sample_id = "consolidated_([0-9]{1,6})_sequence\\.fastq",
  file_sample_id_from_bigwig = "processed_([0-9]{1,6})_sequence_to_S288C_(RPKM|CPM|BPM|RPGC)\\.bw",
  #file_sample_id_from_bigwig = "processed_([0-9]{1,6})_sequence_to_S288C_.*_(RPKM|CPM|BPM|RPGC)\\.bw",
  file_genome_pattern = "S288C_refgenome.fna",
  file_genome_directory = file.path(Sys.getenv("HOME"), "data", "REFGENS"),
  file_feature_directory = file.path(Sys.getenv("HOME"), "data", "feature_files"),
  file_feature_pattern = "eaton_peaks",

  # File Names
  filename_format_group_template = "%s_%s_group%02d_chr%s_%s.svg",
  filename_format_comparison_template = "%s_%s_%s_chr%s_%s.svg",
  title_group_template = paste(
    "%s",         # Title
    "Group: %s",   # Comparison ID
    "Chromosome %s (%d samples)", # Chr info
    "%s",         # Additional info
    "Normalization: %s", # Norm method
    sep = "\n"
  ),
  title_comparison_template = paste(
    "%s",         # Title
    "Comparison: %s",   # Comparison ID
    "Chromosome %s (%d samples)", # Chr info
    "%s",         # Additional info
    "Normalization: %s", # Norm method
    sep = "\n"
  ),
  # Development mode title
  #title_dev_mode = "development",  # Enum: "development" | "publication"
  #title_dev_style = 2,  # Bold
  ## Publication mode title
  #title_pub_template = "%s: Chr%s (%s)",
  #title_pub_size = 1,
  #title_pub_style = 2,  # Bold
  ## Title constraints
  #title_max_width = 40,
  #title_max_lines = 5,
  # Interactive mode
  #interactive_prompt = "Options: [Enter] next plot, 's' skip rest, 'q' quit: ",

  track_defaults_sample = list(
    showaxis = TRUE,
    showtitle = TRUE,
    type = "h",
    size = 1.2,
    background.title = "white",
    fontcolor.title = "black",
    col.border.title = "#e0e0e0",
    cex.title = 0.6,
    fontface = 1,
    title.width = 1.2
  ),

  track_defaults_placeholder = list(
    showaxis = TRUE,
    showtitle = TRUE,
    type = "h",
    size = 0.8,
    background.title = "white",
    background.panel = "#f5f5f5",  # light gray to indicate "empty"
    fontcolor.title = "black",
    col.border.title = "#e0e0e0",
    cex.title = 0.7,
    fontface = 1,
    title.width = 0.9,
    alpha = 0.5,
    grid = FALSE
    #ylim = c(0, 1)          # fixed range for empty tracks
  ),
  track_defaults_control = list(
    showaxis = TRUE,
    showtitle = TRUE,
    type = "h",
    size = 0.8,
    background.title = "white",
    background.panel = "white",
    fontcolor.title = "black",
    col.border.title = "#e0e0e0",
    cex.title = 0.7,
    fontface = 1,
    title.width = 0.9
    #alpha = 0.8
  ),
  track_defaults_feature = list(
    showaxis = FALSE,
    showtitle = TRUE,
    size = 0.5,
    background.title = "white",
    background.panel = "#8b7355",
    fontcolor.title = "black",
    col.border.title = "#e0e0e0",
    cex.title = 0.7,
    fontface = 1,
    title.width = 0.9,
    fill = "#8b4513",
    col = "#8b4513"
  ),

  # New Plot Defaults
  plot_defaults = list(
    margin = 15,
    innerMargin = 5,
    spacing = 10,
    extend.left = 0,
    extend.right = 0,
    col.axis = "black",
    cex.axis = 0.8,
    cex.main = 0.8,
    fontface.main = 2,
    background.panel = "transparent"
  )
)

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
