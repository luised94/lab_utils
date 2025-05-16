#!/usr/bin/env Rscript
# Bootstrap phase
# Also loads OVERRIDE_PRESETS
FUNCTION_FILENAMES <- c("logging", "script_control", "file_operations")
for (function_filename in FUNCTION_FILENAMES) {
    function_filepath <- sprintf("~/lab_utils/core_scripts/functions_for_%s.R", function_filename)
    normalized_path <- normalizePath(function_filepath)
    if (!file.exists(normalized_path)) {
        stop(sprintf("[FATAL] File with functions not found: %s", normalized_path))
    }
    source(normalized_path)
}
message("Bootstrap phase completed...")
if(interactive()) {
  message("Interactive job... sourcing configuration file.")
  script_configuration_path <- "~/lab_utils/core_scripts/script_configuration.R"
  stopifnot(
    "Script configuration file does not exist. Please copy the template." =
    file.exists(script_configuration_path)
  )
  source(script_configuration_path)
  message("Configuration file sourced...")
}

# Ensure the variables expected in the script were //
# defined in the configuration file. //
variables_to_check <- c(
  "EXPERIMENT_IDS",
  "EXPERIMENT_DIR",
  "CHROMOSOMES_TO_PLOT",
  "OUTPUT_FORMAT",
  "ACCEPT_CONFIGURATION",
  "SKIP_PACKAGE_CHECKS"
)
sapply(variables_to_check, function(variable){
  if (!exists(variable)) {
    stop(sprintf(
      fmt = "Variable %s not defined.\nAdjust configuration file.",
      variable
    ))
  }
})
message("All variables defined in the configuration file...")


#} else {
# add the args code?
#}
################################################################################
# Verify Required Libraries
################################################################################
# Add the packages that are used in the script.
required_packages <- c("rtracklayer", "GenomicRanges", "Gviz")
if (!is.character(required_packages) || length(required_packages) == 0) {
  stop("required_packages must be a non-empty character vector")
}

if (!SKIP_PACKAGE_CHECKS) {
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf(
        fmt = "Package '%s' is missing.\nPlease install using renv or base R.",
        pkg
      ))
    }
  }
  SKIP_PACKAGE_CHECKS <- TRUE
}

message("All required packages available...")
################################################################################
# Setup experiment-specific configuration path
################################################################################
number_of_experiments <- length(EXPERIMENT_DIR)
config_paths <- vector("character", length = number_of_experiments)
metadata_paths <- vector("character", length = number_of_experiments)
for (experiment_idx in seq_len(number_of_experiments)) {
  current_experiment_path <- EXPERIMENT_DIR[experiment_idx]
  current_experiment_id <- EXPERIMENT_IDS[experiment_idx]

  config_paths[experiment_idx] <- file.path(
    current_experiment_path, "documentation",
    paste0(current_experiment_id, "_bmc_config.R")
  )
  metadata_paths[experiment_idx] <- file.path(
    current_experiment_path, "documentation",
    paste0(current_experiment_id, "_sample_grid.csv")
  )
}
sapply(c(config_paths, metadata_paths),
  function(file_path){
    if (!file.exists(file_path)) {
      stop(sprintf(
        fmt = "File %s does not exist.\nRun setup_bmc_experiment.R script.",
        file_path
      ))
    }
  }
)

################################################################################
# Setup directories, genome file and file metadata
################################################################################
OUTPUT_DIR <- file.path(EXPERIMENT_DIR[1], "genome_tracks", "final_results")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

metadata_list <- vector("list", length = number_of_experiments)
metadata_categories_list <- vector("list", length = number_of_experiments)

# Patterns for files
bigwig_pattern <- "processed_.*_sequence_to_S288C_blFiltered_CPM\\.bw$"
fastq_pattern <- "consolidated_.*_sequence\\.fastq$"

expected_number_of_samples <- 0
REQUIRED_DIRECTORIES <- c("fastq", "coverage")
for (experiment_idx in seq_len(number_of_experiments)) {
  message("--- Start experiment idx ---")
  current_experiment_path <- EXPERIMENT_DIR[experiment_idx]
  current_config_path <- config_paths[experiment_idx]
  current_metadata_path <- metadata_paths[experiment_idx]
  debug_print(list(
    "title" = "Debug Path information",
    ".Current Experiment path" = current_experiment_path,
    ".Current metadata path" = current_metadata_path,
    ".Current Config path" = current_config_path
  ))

  # Build all required directory paths for this experiment
  required_paths <- file.path(current_experiment_path, REQUIRED_DIRECTORIES)
  names(required_paths) <- REQUIRED_DIRECTORIES
  missing_dirs <- required_paths[!dir.exists(required_paths)]
  if (length(missing_dirs) > 0) {
    stop("Missing required experiment subdirectories: ", 
         paste(missing_dirs, collapse = ", "))
  }

  # Load *_CONFIG variables
  source(current_config_path)
  current_metadata_df <- read.csv(current_metadata_path, stringsAsFactors = FALSE)

  current_metadata_df$experiment_id <- EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID
  metadata_categories_list[[experiment_idx]] <- EXPERIMENT_CONFIG$CATEGORIES
  #current_metadata_df$bigwig_file_paths <- EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID
  # Add determination of sample ids and addition to metadata frame
  metadata_list[[experiment_idx]] <- current_metadata_df
  expected_number_of_samples <- expected_number_of_samples + EXPERIMENT_CONFIG$METADATA$EXPECTED_SAMPLES
  message("--- End iteration ---")
  message("\n")
} # end for loop

# Get the union of the categories
all_keys <- unique(unlist(lapply(metadata_categories_list, names)))
merged_categories <- setNames(
  lapply(all_keys, function(key) {
    unique(unlist(lapply(metadata_categories_list, function(lst) lst[[key]])))
  }),
  all_keys
)

# Get union of all column names
all_cols <- unique(unlist(lapply(metadata_list, names)))
# Reorder and fill missing columns for each dataframe
metadata_aligned <- lapply(
  metadata_list,
  function(df) {
    # Add missing cols as NA
    missing <- setdiff(all_cols, names(df))
    df[missing] <- NA
    # Reorder columns
    df <- df[all_cols]
    return(df)
  }
)
# Now bind
metadata_df <- do.call(rbind, metadata_aligned)

#####################
# Load the metadata dataframe with all experiments.
#####################
