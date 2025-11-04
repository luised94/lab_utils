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
ROOT_DIRECTORY <- system("git rev-parse --show-toplevel", intern = TRUE)
if (length(ROOT_DIRECTORY) == 0 || ROOT_DIRECTORY == "") {
  stop("Not in git repository. Current directory: ", getwd())
}
CORE_SCRIPTS_PATH <- file.path(ROOT_DIRECTORY, "core_scripts")
SETUP_ENV_SCRIPT_PATH <- file.path(CORE_SCRIPTS_PATH, "setup_ngs_script_environment.R")
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
  "CHROMOSOMES_TO_PLOT",
  "OUTPUT_EXTENSION", "BIGWIG_PATTERN",
  "FASTQ_PATTERN", "SAMPLE_ID_PATTERN",
)

is_missing_variable <- !sapply(required_configuration_variables, exists)
missing_variables <- required_configuration_variables[is_missing_variable]

if (length(missing_vars) > 0) {
  stop(
    "Missing required configuration variables in 'configuration_script_bmc.R':\n",
    "  ", paste(missing_vars, collapse = "\n  ")
  )
}

message("All variables defined in the configuration file...")

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
