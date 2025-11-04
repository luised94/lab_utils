################################################################################
# Prepare enriched metadata csv to use as entry point to other scripts
# Author: Luis | Date: 2024-11-27 | Version: 2.0.0
################################################################################
# PURPOSE:
#   Creates standardized directory structure and metadata files for BMC ChIP-seq experiments
#
# USAGE:
#   1) Update configuration_experiment_bmc.R file.
#   2) From R REPL: source("core_scripts/plot_genome_tracks_for_multiple_chromosomes_experiments_an.R")

#
# DEPENDENCIES:
#   configuration_experiment_bmc.R
#   configuration_script_bmc.R
#   ~/lab_utils/core_scripts/setup_bmc_experiment.R outputs
#   ~/lab_utils/core_scripts/{logging,script_control,file_operations}.R
#   required_packages
#
# OUTPUTS:
# - Svg or pdf files with genome tracks from multiple experiment ids for chromosomes and categories.
################################################################################
#===============================================================================
# BOOTSTRAP:
#   1) Determine Repository Root
#   2) Initialize environment via script.
#   3) Includes sourcing configuration_script_bmc.R.
#===============================================================================
message("=== BOOTSTRAP: Project root and environment ===")

# Use git to set root directory.
# Ensure different worktrees work from different directories.
ROOT_DIRECTORY <- system("git rev-parse --show-toplevel", intern = TRUE)
if (length(ROOT_DIRECTORY) == 0 || ROOT_DIRECTORY == "") {
  stop("Not in git repository. Current directory: ", getwd())
}
ROOT_DIRECTORY <- normalizePath(ROOT_DIRECTORY, mustWork = TRUE)

CORE_SCRIPTS_PATH <- file.path(ROOT_DIRECTORY, "core_scripts")
SETUP_ENV_SCRIPT_PATH <- file.path(CORE_SCRIPTS_PATH, "setup_ngs_script_environment.R")
if (!dir.exists(CORE_SCRIPTS_PATH)) {
  stop("Missing directory with scripts: ", SETUP_ENV_SCRIPT_PATH)
}

if (!file.exists(SETUP_ENV_SCRIPT_PATH)) {
  stop("Environment setup script missing: ", SETUP_ENV_SCRIPT_PATH)
}

message("Repository root path: ", ROOT_DIRECTORY)
message("Setup env script path: ", SETUP_ENV_SCRIPT_PATH)
message("Sourcing environment setup script...")
source(SETUP_ENV_SCRIPT_PATH)
message("Reentering analysis script...")

# Ensure the variables expected in the script were //
# defined in the configuration file. //
# configuration_script_bmc.R //
required_configuration_variables <- c(
  "EXPERIMENT_IDS", "EXPERIMENT_DIR_PATHS",
  "CHROMOSOMES_TO_PLOT", "OUTPUT_EXTENSION",
  "BIGWIG_PATTERN", "FASTQ_PATTERN", "SAMPLE_ID_PATTERN",
)

is_missing_variable <- !sapply(required_configuration_variables, exists)
missing_variables <- required_configuration_variables[is_missing_variable]

if (length(missing_vars) > 0) {
  stop(
    "Missing required configuration variables in 'configuration_script_bmc.R':\n",
    "  ", paste(missing_vars, collapse = "\n  ")
  )
}

# Enforce single-experiment mode for this script
# To run for multiple experiments, in R console:
#for (exp_id in c("241105Bel", "241106Bel", "241107Bel")) {
#  EXPERIMENT_IDS <- exp_id
#  source("core_scripts/prepare_experiment_metadata.R")
#}
if (length(EXPERIMENT_IDS) != 1) {
  stop(
    "This script processes one experiment at a time.\n",
    "Set EXPERIMENT_IDS to a single ID in configuration.\n",
    "Found: ", paste(EXPERIMENT_IDS, collapse = ", ")
  )
}

EXPERIMENT_ID <- EXPERIMENT_IDS[1]  # Simplify for single-experiment context

message("Bootstrap complete. Processing experiment: ", EXPERIMENT_ID, "\n")

#===============================================================================
# SECTION 1: Load Required Packages
#===============================================================================
message("=== SECTION 1: Load Required Packages ===")

required_packages <- c("rtracklayer", "GenomicRanges", "Biostrings")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package not available: ", pkg, "\n",
         "Install with: BiocManager::install('", pkg, "')")
  }
  library(pkg, character.only = TRUE)
}

message("Packages loaded: ", paste(required_packages, collapse = ", "))
message("Package verification complete...")


#===============================================================================
# SECTION 2: Select and Validate Experiment
#===============================================================================
message("=== SECTION 2: Select and Validate Experiment ===")

# Handle multiple EXPERIMENT_IDS: process first only
if (length(EXPERIMENT_IDS) > 1) {
  message("Multiple experiment IDs detected: ", paste(EXPERIMENT_IDS, collapse = ", "))
  message("Processing first experiment only: ", EXPERIMENT_IDS[1])
  EXPERIMENT_ID <- EXPERIMENT_IDS[1]
  EXPERIMENT_DIR_PATH <- EXPERIMENT_DIR_PATHS[1]
} else {
  EXPERIMENT_ID <- EXPERIMENT_IDS[1]
  EXPERIMENT_DIR_PATH <- EXPERIMENT_DIR_PATHS[1]
}

# Validate experiment directory exists
if (!dir.exists(EXPERIMENT_DIR_PATH)) {
  stop("Experiment directory not found: ", EXPERIMENT_DIR_PATH, "\n",
       "Run setup_bmc_experiment.R first for experiment: ", EXPERIMENT_ID)
}

# Build paths to required files
experiment_config_path <- file.path(
  EXPERIMENT_DIR_PATH,
  "documentation",
  paste0(EXPERIMENT_ID, "_configuration_experiment_bmc.R")
)

sample_grid_path <- file.path(
  EXPERIMENT_DIR_PATH,
  "documentation",
  paste0(EXPERIMENT_ID, "_sample_grid.csv")
)

# Validate required files exist
if (!file.exists(experiment_config_path)) {
  stop("Experiment configuration not found: ", basename(experiment_config_path), "\n",
       "Expected location: ", experiment_config_path)
}

if (!file.exists(sample_grid_path)) {
  stop("Sample grid not found: ", basename(sample_grid_path), "\n",
       "Run setup_bmc_experiment.R first")
}

message("Experiment ID: ", EXPERIMENT_ID)
message("Experiment directory: ", EXPERIMENT_DIR_PATH)
message("Configuration: ", basename(experiment_config_path))
message("Sample grid: ", basename(sample_grid_path))
message("Section 2 complete\n")
#===============================================================================
# SECTION 3: Load Experiment Configuration
#===============================================================================
message("=== SECTION 3: Load Experiment Configuration ===")

# Source experiment-specific configuration
source(experiment_config_path)

# Validate EXPERIMENT_CONFIG exists
if (!exists("EXPERIMENT_CONFIG") || !is.list(EXPERIMENT_CONFIG)) {
  stop("EXPERIMENT_CONFIG not defined or not a list.\n",
       "Check file: ", basename(experiment_config_path))
}

# Validate EXPERIMENT_CONFIG has required grouping columns
#if (!exists("EXPERIMENTAL_GROUPING_COLUMNS", where = EXPERIMENT_CONFIG)) {
#  stop("EXPERIMENT_CONFIG missing EXPERIMENTAL_GROUPING_COLUMNS.\n",
#       "Check configuration file")
#}

message("Experiment configuration loaded")
message("Section 3 complete\n")
#===============================================================================
# SECTION 4: Load Sample Grid
#===============================================================================
message("=== SECTION 4: Load Sample Grid ===")

# Load sample grid CSV
metadata_df <- read.csv(sample_grid_path, stringsAsFactors = FALSE)

# Basic validation
if (nrow(metadata_df) == 0) {
  stop("Loaded metadata grid has zero rows: ", basename(sample_grid_path))
}

message("Sample grid loaded")
message("  Samples: ", nrow(metadata_df))
message("  Columns: ", ncol(metadata_df))
message("  Column names: ", paste(colnames(metadata_df), collapse = ", "))
message("Section 4 complete\n")
