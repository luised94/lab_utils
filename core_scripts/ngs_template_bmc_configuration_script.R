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

#===============================================================================
# FILE PATTERNS
#===============================================================================
FASTQ_PATTERN <- "fastpfiltered_.*_(R1|R2|NA)\\.fastq$"
SAMPLE_ID_PATTERN <- "fastpfiltered_(D[0-9]{2}-[0-9]{1,6})_(R1|R2|NA)\\.fastq$"
EXPERIMENT_ID_PATTERN <- "^\\d{6}Bel$"
GENOME_FILE_PATTERN <- "S288C_refgenome.fna"
FEATURE_FILE_PATTERN <- "eaton_peaks"

#===============================================================================
# OUTPUT SETTINGS
#===============================================================================
OUTPUT_FORMAT <- "svg"
OUTPUT_FORMATS_VALID <- c("svg", "pdf", "png")

#===============================================================================
# PATHS
#===============================================================================
GENOME_DIR_PATH <- file.path(Sys.getenv("HOME"), "data", "REFGENS")
FEATURE_DIR_PATH <- file.path(Sys.getenv("HOME"), "data", "feature_files")


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

#===============================================================================
# GENOME TRACK VISUALIZATION
#===============================================================================
CHROMOSOMES_TO_PLOT <- c(7, 10, 14)
BIGWIG_NORM_METHOD <- "CPM"

# Display dimensions
GTRACK_DISPLAY_WIDTH <- 10
GTRACK_DISPLAY_HEIGHT <- 8

# Track parameters
GTRACK_YLIM_DEFAULT <- c(0, 1000)
GTRACK_Y_AXIS_SCALING <- "individual"  # Changed from array to single value
#TRACK_SCALING_MODES_VALID <- c("local", "individual")

# Track colors
GTRACK_COLOR_INPUT <- "#808080"

# Track naming formats
GTRACK_NAME_FORMAT_SAMPLE <- "%s: %s"
GTRACK_NAME_FORMAT_CONTROL <- "%s: %s - %s"
GTRACK_NAME_FORMAT_PLACEHOLDER <- "%s: %s - %s"
GTRACK_NAME_FORMAT_SUFFIX <- "(No data)"
GTRACK_NAME_FORMAT_AXIS <- "Chr %s Axis"

# Filename formats
GTRACK_FILENAME_GROUP <- "%s_%s_group%02d_chr%s_%s.svg"
GTRACK_FILENAME_COMPARISON <- "%s_%s_%s_chr%s_%s.svg"

# Title templates (multi-line)
GTRACK_TITLE_GROUP <- paste(
  "%s",
  "Group: %s",
  "Chromosome %s (%d samples)",
  "%s",
  "Normalization: %s",
  sep = "\n"
)

GTRACK_TITLE_COMPARISON <- paste(
  "%s",
  "Comparison: %s",
  "Chromosome %s (%d samples)",
  "%s",
  "Normalization: %s",
  sep = "\n"
)

# Styling parameters for different track types

GTRACK_DEFAULTS_SAMPLE <- list(
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
)

GTRACK_DEFAULTS_FEATURE <- list(
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
)

GTRACK_PLOT_DEFAULTS <- list(
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
