################################################################################
# Title: Setup NGS Script Environment
# Author: Luis/Claude/Qwen | Date: 2025-11-03 | Version: 1.0
################################################################################
# Purpose:
#   Initialize analysis environment by sourcing configuration, processing
#   experiment IDs, computing runtime values, and validating base paths.
#   Called by entry point scripts after bootstrap block.
#
# Prerequisites (must be defined before sourcing):
#   ROOT_DIRECTORY      - Git repository root (from git rev-parse)
#   CORE_SCRIPTS_PATH   - Path to core_scripts/ directory
#
# Usage:
#   # In entry point script (see bootstrap block):
#   ROOT_DIRECTORY <- system("git rev-parse --show-toplevel", intern = TRUE)
#   CORE_SCRIPTS_PATH <- file.path(ROOT_DIRECTORY, "core_scripts")
#   source(file.path(CORE_SCRIPTS_PATH, "setup_ngs_script_environment.R"))
#
# Dependencies:
#   - SOURCES: configuration_script_bmc.R (in CORE_SCRIPTS_PATH)
#   - ACCESSES: HOME environment variable
#
# Output:
#   Global variables:
#     EXPERIMENT_IDS        - Character vector of experiment IDs (length >= 1)
#     EXPERIMENT_DIR_PATHS        - Named vector of experiment directory paths
#     CURRENT_TIMESTAMP     - Runtime timestamp (YYYYMMDD_HHMMSS)
#     CURRENT_DATE          - Runtime date (YYYYMMDD)
#     OUTPUT_EXTENSION      - File extension for plots (e.g., ".svg")
#     BIGWIG_PATTERN        - Regex pattern for bigwig files
#
# Common Failure Points:
#   1. ROOT_DIRECTORY not defined -> calling script must set it first
#   2. configuration_script_bmc.R missing -> verify file exists
#   3. HOME environment variable not set -> check system environment
#   4. Invalid EXPERIMENT_IDS format -> check EXPERIMENT_ID_PATTERN
#   5. GENOME_DIR_PATH or FEATURE_DIR_PATH missing -> verify paths in config
################################################################################
#===============================================================================
# SECTION 1: Validate Prerequisites
#===============================================================================
message("Environment setup: Validating prerequisites...")

# Verify ROOT_DIRECTORY was set by calling script
if (!exists("ROOT_DIRECTORY") || is.null(ROOT_DIRECTORY) || ROOT_DIRECTORY == "") {
  stop("ROOT_DIRECTORY not defined.\n",
       "Calling script must set ROOT_DIRECTORY before sourcing this script")
}

# Verify CORE_SCRIPTS_PATH was set by calling script
if (!exists("CORE_SCRIPTS_PATH") || !dir.exists(CORE_SCRIPTS_PATH)) {
  stop("CORE_SCRIPTS_PATH not defined or does not exist: ", CORE_SCRIPTS_PATH)
}

message("Prerequisites validated")

#===============================================================================
# SECTION 2: Source Configuration Files
#===============================================================================
message("Environment setup: Sourcing configuration...")

# Build path to script configuration
SCRIPT_CONFIG_PATH <- file.path(CORE_SCRIPTS_PATH, "configuration_script_bmc.R")

# Validate configuration file exists
if (!file.exists(SCRIPT_CONFIG_PATH)) {
  stop("Script configuration file not found: ", SCRIPT_CONFIG_PATH, "\n",
       "Check configuration_script_bmc.R exists in core_scripts/")
}

# Source configuration
source(SCRIPT_CONFIG_PATH)

# Validate required environment variables
home_dir <- Sys.getenv("HOME")
if (home_dir == "") {
  stop("HOME environment variable not set.\n",
       "Check system environment configuration")
}

# Validate sentinel configuration variables loaded
required_config_vars <- c(
  "EXPERIMENT_IDS",
  "GENOME_DIR_PATH",
  "FEATURE_DIR_PATH",
  "OUTPUT_FORMAT",
  "BIGWIG_NORM_METHOD",
  "EXPERIMENT_ID_PATTERN",
  "TIMESTAMP_FORMAT",
  "DATE_FORMAT"
)

missing_vars <- required_config_vars[!sapply(required_config_vars, exists)]

if (length(missing_vars) > 0) {
  stop("Missing required configuration variables:\n",
       paste(missing_vars, collapse = ", "), "\n",
       "Check configuration_script_bmc.R")
}

message("Configuration loaded")

#===============================================================================
# SECTION 3: Process EXPERIMENT_IDS
#===============================================================================
message("Environment setup: Processing experiment IDs...")

# Check if EXPERIMENT_IDS contains comma-separated values
is_comma_separated <- grepl(",", EXPERIMENT_IDS)

if (is_comma_separated) {
  # Split and clean
  experiment_ids_vec <- strsplit(EXPERIMENT_IDS, ",", fixed = TRUE)[[1]]
  experiment_ids_vec <- trimws(experiment_ids_vec)
  experiment_ids_vec <- experiment_ids_vec[experiment_ids_vec != ""]

  # Validate each ID format
  for (exp_id in experiment_ids_vec) {
    if (!grepl(EXPERIMENT_ID_PATTERN, exp_id, perl = TRUE)) {
      stop("Invalid experiment ID format: ", exp_id, "\n",
           "Expected pattern: ", EXPERIMENT_ID_PATTERN)
    }
  }

  # Check for duplicates
  if (any(duplicated(experiment_ids_vec))) {
    duplicate_ids <- experiment_ids_vec[duplicated(experiment_ids_vec)]
    stop("Duplicate experiment IDs detected: ",
         paste(duplicate_ids, collapse = ", "))
  }

  # Update EXPERIMENT_IDS with cleaned vector
  EXPERIMENT_IDS <- experiment_ids_vec

} else {
  # Validate single ID format
  if (!grepl(EXPERIMENT_ID_PATTERN, EXPERIMENT_IDS, perl = TRUE)) {
    stop("Invalid experiment ID format: ", EXPERIMENT_IDS, "\n",
         "Expected pattern: ", EXPERIMENT_ID_PATTERN)
  }

  # Convert to vector of length 1 for consistent handling
  EXPERIMENT_IDS <- c(EXPERIMENT_IDS)
}

# Create parallel array of experiment directory paths
EXPERIMENT_DIR_PATHS <- file.path(home_dir, "data", EXPERIMENT_IDS)
names(EXPERIMENT_DIR_PATHS) <- EXPERIMENT_IDS

# Check which directories exist
existing_dirs <- EXPERIMENT_DIR_PATHS[dir.exists(EXPERIMENT_DIR_PATHS)]
missing_dirs <- EXPERIMENT_DIR_PATHS[!dir.exists(EXPERIMENT_DIR_PATHS)]

if (length(missing_dirs) > 0) {
  warning("Some experiment directories do not exist:\n",
          paste(missing_dirs, collapse = "\n"), "\n",
          "These may need to be created with setup_bmc_experiment.R")
}

message("Experiment IDs: ", paste(EXPERIMENT_IDS, collapse = ", "))
message("Existing directories: ", length(existing_dirs), " of ", length(EXPERIMENT_IDS))

#===============================================================================
# SECTION 4: Compute Runtime Values
#===============================================================================
message("Environment setup: Computing runtime values...")

# Compute timestamp and date using formats from config
CURRENT_TIMESTAMP <- format(Sys.time(), TIMESTAMP_FORMAT)
CURRENT_DATE <- format(Sys.Date(), DATE_FORMAT)

# Compute derived constants from configuration
OUTPUT_EXTENSION <- paste0(".", OUTPUT_FORMAT)
BIGWIG_PATTERN <- sprintf("D[0-9]{2}-[0-9]{1,6}_%s\\.bw", BIGWIG_NORM_METHOD)

message("Timestamp: ", CURRENT_TIMESTAMP)
message("Output format: ", OUTPUT_FORMAT)

#===============================================================================
# SECTION 5: Validate Base Paths
#===============================================================================
message("Environment setup: Validating base paths...")

# Validate genome directory exists
if (!dir.exists(GENOME_DIR_PATH)) {
  stop("GENOME_DIR_PATH not found: ", GENOME_DIR_PATH, "\n",
       "Check GENOME_DIR_PATH in configuration_script_bmc.R")
}

# Validate feature directory exists
if (!dir.exists(FEATURE_DIR_PATH)) {
  stop("FEATURE_DIR_PATH not found: ", FEATURE_DIR_PATH, "\n",
       "Check FEATURE_DIR_PATH in configuration_script_bmc.R")
}

message("Genome directory: ", GENOME_DIR_PATH)
message("Feature directory: ", FEATURE_DIR_PATH)
message("Environment setup complete")

################################################################################
# End of setup_ngs_script_environment.R
################################################################################
