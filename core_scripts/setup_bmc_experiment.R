#!/usr/bin/env Rscript
################################################################################
# BMC ChIP-seq Experiment Setup
# Author: Luis | Date: 2025-10-25 | Version: 2.1.0
################################################################################
# PURPOSE:
#   Creates standardized directory structure and metadata files for BMC ChIP-seq experiments
#
# USAGE:
#   1) Update configuration_experiment_bmc.R file.
#   2) From R REPL: source("core_scripts/setup_bmc_experiment.R")
#
# DEPENDENCIES: ~/lab_utils/core_scripts/configuration_experiment_bmc.R
#
# OUTPUTS:
# - Standard directory structure (peak/, fastq/, alignment/, bigwig/, plots/)
# - Sample metadata files (_sample_grid.csv, _bmc_table.tsv)
################################################################################
if(interactive()) {
  message("Running from repl... Loading functions.")
} else {
  stop("Run the script from the R repl in an interactive session.")
}

ROOT_DIRECTORY <- system("git rev-parse --show-toplevel", intern = TRUE)
CORE_SCRIPTS_PATH <- file.path(ROOT_DIRECTORY, "core_scripts")
SCRIPT_CONFIGURATION_PATH <- file.path(CORE_SCRIPTS_PATH, "configuration_script_bmc.R")
EXPERIMENT_CONFIGURATION_PATH <- file.path(
  CORE_SCRIPTS_PATH,
  "configuration_experiment_bmc.R"
)
FUNCTION_FILENAME_TEMPLATE <- file.path(CORE_SCRIPTS_PATH, "functions_for_%s.R")

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

#---------------------------------------
# Load user functions
#---------------------------------------
function_filenames <- c(
  "logging", "script_control",
  "file_operations", "bmc_config_validation"
)

for (function_filename in function_filenames) {
  function_filepath <- sprintf(
    FUNCTION_FILENAME_TEMPLATE, 
    function_filename
  )
  normalized_path <- normalizePath(function_filepath)
  if (!file.exists(normalized_path)) {
    stop(sprintf("[FATAL] File with functions not found: %s", normalized_path))
  }
  source(normalized_path)
}

message("Loaded functions... Sourcing configuration.")

# Ensure configuration files exist.
stopifnot(
  "Script configuration file does not exist. Please copy the template." =
    file.exists(SCRIPT_CONFIGURATION_PATH),
  "Experiment configuration file does not exist. Please copy the template." =
    file.exists(EXPERIMENT_CONFIGURATION_PATH)
)

#-------------------------------------------------------------------------------
# Validate experiment configuration
#-------------------------------------------------------------------------------
source(EXPERIMENT_CONFIGURATION_PATH)

stopifnot(
  "EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID does not match configuration EXPERIMENT_ID" =
    identical(EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID, EXPERIMENT_IDS),
  "Experiment configuration file does not exist. Please copy the template." =
    file.exists(EXPERIMENT_CONFIGURATION_PATH)
)

EXPERIMENT_IDS <- EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID
category_names <- names(EXPERIMENT_CONFIG$CATEGORIES)

stopifnot(
  "Experiment ID must be a character string" =
    is.character(EXPERIMENT_IDS),
  "Invalid experiment ID format. Expected: YYMMDD'Bel'" =
    grepl("^\\d{6}Bel$", EXPERIMENT_IDS)
)

required_sections <- c(
  "METADATA", "CATEGORIES", "INVALID_COMBINATIONS",
  "COMPARISONS", "CONTROL_FACTORS",
  "COLUMN_ORDER", "NORMALIZATION"
  #"EXPERIMENTAL_CONDITIONS","SAMPLE_CLASSIFICATIONS"
)

missing_sections <- setdiff(required_sections, names(EXPERIMENT_CONFIG))
if (length(missing_sections) > 0) {
  stop(sprintf("Missing required config sections: %s",
        paste(missing_sections, collapse = ", ")))
}

# Validate configuration structure
stopifnot(
  "Missing required config sections" =
    all(required_sections %in% names(EXPERIMENT_CONFIG))
)

sapply(category_names,
  function(category_name) {
    category_values <- EXPERIMENT_CONFIG$CATEGORIES[[category_name]]
    is_not_character <- !is.character(category_values)
    is_duplicated <- duplicated(category_values)
    if (any(is_not_character)) {
      stop(sprintf("Some category values in '%s' category are not character: %s",
        category_name,
        paste(category_values[is_not_character], collapse = ", ")
      ))
    }
    if (any(is_duplicated)) {
      stop(sprintf("Some category values in '%s' category are duplicated: %s",
        category_name,
        paste(category_values[is_duplicated], collapse = ", ")
      ))
    }
  }
)

categories_and_column_order_are_not_identical <- !identical(
  sort(category_names),
  sort(EXPERIMENT_CONFIG$COLUMN_ORDER)
)

if (categories_and_column_order_are_not_identical) {
  stop("Column order must include all category columns.")
}

lapply(names(EXPERIMENT_CONFIG$INVALID_COMBINATIONS),
  function(invalid_combination_condition) {
    combination_condition <- EXPERIMENT_CONFIG$INVALID_COMBINATIONS[[invalid_combination_condition]]
    invalid_category <- setdiff(
      all.vars(combination_condition),
      category_names
    )
    if (length(invalid_category) > 0) {
      stop(sprintf(
        "Invalid columns in INVALID_COMBINATIONS '%s':\n%s",
        combination_condition, paste(invalid_category, collapse = ", ")
      ))

    }
  }
)

lapply(names(EXPERIMENT_CONFIG$CONTROL_FACTORS),
  function(control_factor_names) {
    control_factor_categories <- EXPERIMENT_CONFIG$CONTROL_FACTORS[[control_factor_names]]
    invalid_category <- setdiff(
      control_factor_categories,
      category_names
    )
    if (length(invalid_category) > 0) {
      stop(sprintf(
        "Invalid columns in CONTROL_FACTORS '%s':\n%s",
        control_factor_names, paste(invalid_category, collapse = ", ")
      ))

    }
  }
)

cat("\n[VALIDATED] Experiment configuration loaded successfully\n")

source(SCRIPT_CONFIGURATION_PATH)
message("Configuration file sourced... Checking configuration variables.")

# Ensure the variables expected in the script were //
# defined in the configuration file. //
# configuration_script_bmc.R //
required_configuration_variables <- c(
  "EXPERIMENT_IDS", "EXPERIMENT_DIR",
  "ACCEPT_CONFIGURATION"
)

is_missing_variable <- !sapply(required_configuration_variables, exists)
missing_variables <- required_configuration_variables[is_missing_variable]

if (length(missing_variables) > 0 ) {
  stop("Missing variable. Please define in 'configuration_script_bmc.R' file.",
     paste(missing_variables, collapse = ", "))
}
stopifnot(
  "Only one experiment id required for this script" = length(EXPERIMENT_IDS) == 1
)

message("All variables defined in the configuration file...")

required_configs <- c("EXPERIMENT_CONFIG", "RUNTIME_CONFIG")
validate_configs(required_configs)
invisible(lapply(required_configs, function(config) {
  print_config_settings(get(config), title = config)
}))

stopifnot(
  "Script EXPERIMENT_IDS is not the same as CONFIG EXPERIMENT_IDS" =
    EXPERIMENT_IDS == EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID
)

message("Experiment configuration...")

#-------------------------------------------------------------------------------
# Directory Setup
#-------------------------------------------------------------------------------
# Create directory structure
full_paths <- file.path(EXPERIMENT_DIR, DATA_DIRECTORIES)
invisible(lapply(full_paths, function(path) {
  if (RUNTIME_CONFIG$output_dry_run) {
    cat(sprintf("[DRY RUN] Would create directory: %s\n", path))
  } else {
    dir_created <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (RUNTIME_CONFIG$debug_verbose) {
      status <- if (dir_created) "Created" else "Already exists"
      cat(sprintf("[%s] %s\n", status, path))
    }
  }
}))

# Report directory creation status
if (RUNTIME_CONFIG$debug_verbose) {
  mode <- if (RUNTIME_CONFIG$output_dry_run) "DRY RUN" else "LIVE RUN"
  cat(sprintf("\n[%s] Directory structure for experiment: %s\n", mode, EXPERIMENT_IDS))
  cat(sprintf("[%s] Base directory: %s\n", mode, EXPERIMENT_DIR))
}

cat("Directories created successfully!\n")
#-------------------------------------------------------------------------------
# Metadata Generation, Organization and Validation
#-------------------------------------------------------------------------------
# Generate experimental combinations
metadata <- do.call(expand.grid, EXPERIMENT_CONFIG$CATEGORIES)

# Filter invalid combinations if set.
# Run without combination to see outputs.
if (length(EXPERIMENT_CONFIG$INVALID_COMBINATIONS) > 0) {
  invalid_idx <- Reduce(
    `|`,
    lapply(EXPERIMENT_CONFIG$INVALID_COMBINATIONS, eval, envir = metadata)
  )
  metadata <- subset(metadata, !invalid_idx)
}

# Apply experimental conditions
#valid_idx <- Reduce(
#  `|`,
#  lapply(EXPERIMENT_CONFIG$EXPERIMENTAL_CONDITIONS, eval, envir = metadata)
#)
#metadata <- subset(metadata, valid_idx)

# Enforce factor levels from config
for (col_name in names(EXPERIMENT_CONFIG$CATEGORIES)) {
  if (col_name %in% colnames(metadata)) {
    metadata[[col_name]] <- factor(
      metadata[[col_name]],
      levels = EXPERIMENT_CONFIG$CATEGORIES[[col_name]],
      ordered = TRUE
    )
  }
}

# Sort metadata according to column order
metadata <- metadata[
  do.call(
    order,
    metadata[EXPERIMENT_CONFIG$COLUMN_ORDER]
  ),
]

# Verify sample count
n_samples <- nrow(metadata)
expected <- EXPERIMENT_CONFIG$METADATA$EXPECTED_SAMPLES
if (n_samples != expected) {
  terminal_width <- get_width()
  # Set width to 200
  options(width = terminal_width)
  # Print diagnostic information
  cat("\nDiagnostic Information:\n")
  cat("----------------------\n")
  print(metadata)
  stop(sprintf("Expected %d samples, got %d", expected, n_samples))
}

#-------------------------------------------------------------------------------
# Adjust order of metadata
#-------------------------------------------------------------------------------
# If mistake was made during submission, this section of code will readjust based on the ROW_ORDER_CORRECTION value
#if(is.null(ROW_ORDER_CORRECTION)) {
#
#}
  #print(table(metadata$antibody))  # Show antibody distribution
  #cat("\nFull sample breakdown:\n")
  #print(summary(metadata))     # Show all category distributions
  #cat("\n")
  ## Control this in the configuration_experiment_bmc.R file.
  ## Helps display metadata values to help narrow down where you need to add combinations to filter.
  #if (show_all_metadata) {
  #  print(metadata)
  #} else if (show_particular_metadata) {
  #  print(metadata[metadata[, category_to_show] == values_to_show, ])
  #}

#-------------------------------------------------------------------------------
# Sample Classification
#-------------------------------------------------------------------------------
#sample_classifications <- EXPERIMENT_CONFIG$SAMPLE_CLASSIFICATIONS
#
## First, create a matrix/data frame to store all classification results
#classification_results <- matrix(FALSE,
#                 nrow = nrow(metadata),
#                 ncol = length(sample_classifications),
#                 dimnames = list(NULL, names(sample_classifications)))
#
## Evaluate each classification condition
#for (type in names(sample_classifications)) {
#  classification_results[, type] <- eval(sample_classifications[[type]],
#                     envir = metadata)
#}
#
## Create the final classification vector
## Default classification
#metadata$sample_type <- "treatment"
#for (type in names(sample_classifications)) {
#  # Find rows where this classification is TRUE
#  matching_rows <- classification_results[, type]
#  # Assign the type name (removing 'is_' prefix)
#  metadata$sample_type[matching_rows] <- sub("^is_", "", type)
#}
#
## Validation check
#multiple_classifications <- rowSums(classification_results) > 1
#if (any(multiple_classifications)) {
#  cat("\nERROR: Multiple Classification Detected!\n")
#  cat("----------------------------------------\n")
#
#  # Show problematic samples with their classifications
#  problem_samples <- metadata[multiple_classifications, ]
#  cat("Samples with multiple classifications:\n\n")
#
#  # Show which classifications were TRUE for each problematic sample
#  for (i in which(multiple_classifications)) {
#    cat(sprintf("\nSample %d:\n", i))
#    cat("Sample details:\n")
#    print(metadata[i, ])
#    cat("\nMatching classifications:\n")
#    matching_types <- names(classification_results[i,])[classification_results[i,]]
#    print(matching_types)
#    cat("----------------------------------------\n")
#  }
#
#  stop("Please fix multiple classifications in experiment configuration")
#}
#
## Success diagnostic display
#cat("\nSample Classification Summary:\n")
#cat("============================\n")
#
## Overall counts
#cat("\n1. Distribution of sample types:\n")
#print(table(metadata$sample_type))
#
## Detailed breakdown by relevant factors
#cat("\n2. Sample types by antibody:\n")
#print(table(metadata$sample_type, metadata$antibody))
#
## Show a few samples from each classification
#cat("\n3. Example samples from each classification:\n")
#for (type in unique(metadata$sample_type)) {
#  cat(sprintf("\n%s samples:\n", toupper(type)))
#  print(metadata[metadata$sample_type == type, ][1:min(3, sum(metadata$sample_type == type)), ])
#  cat("----------------------------------------\n")
#}
#
## Verification message
#cat("\nClassification Verification:\n")
#cat(sprintf("- Total samples: %d\n", nrow(metadata)))
#cat(sprintf("- Classified samples: %d\n", sum(table(metadata$sample_type))))
#cat(sprintf("- Unclassified samples: %d\n", sum(is.na(metadata$sample_type))))

# Generate sample names
metadata$full_name <- apply(metadata, 1, paste, collapse = "_")
metadata$short_name <- apply(
  metadata[, EXPERIMENT_CONFIG$COLUMN_ORDER],
  1,
  function(x) paste0(substr(x, 1, 1), collapse = "")
)

#-------------------------------------------------------------------------------
# BMC Metadata Generation
#-------------------------------------------------------------------------------
bmc_metadata <- data.frame(
  SampleName = metadata$full_name,
  Vol_uL = 15,
  Conc = 0,
  Type = ifelse(
    tolower(metadata$antibody) == "input",
    "Input",
    "ChIP"
  ),
  Genome = "Saccharomyces cerevisiae",
  Notes = ifelse(
    tolower(metadata$antibody) == "input",
    "Run on fragment analyzer. Usually diluted 2 fold.",
    "Run on femto pulse. Samples usually run without dilution."
  ),
  Pool = "A",
  stringsAsFactors = FALSE
)

#-------------------------------------------------------------------------------
# File Output Generation
#-------------------------------------------------------------------------------
if (ACCEPT_CONFIGURATION) {
  "Accept configuration enabled... Continuing"
} else {
  stop("Configuration not accepted... Update ACCEPT_CONFIGURATION when ready to continue")
}
filenames <- c("sample_grid.csv", "bmc_table.tsv", "configuration_experiment_bmc.R")
# Loop through each filename to handle path assignment and file writing
for (filename in filenames) {
  # Construct the output file path
  output_file_path <- file.path(
    EXPERIMENT_DIR,
    "documentation",
    paste0(EXPERIMENT_IDS, "_", filename)
  )
  # Handle file writing with dry run checks
  if (RUNTIME_CONFIG$output_dry_run) {
    # Dry-run message
    if (file.exists(output_file_path)) {
      cat(sprintf("[DRY RUN] File exists: %s. Overwrite will occur if dry-run is disabled.\n", output_file_path))
    } else {
      cat(sprintf("[DRY RUN] Would write file to: %s\n", output_file_path))
    }
  } else {
    # Determine the appropriate write function based on the file extension
    if (endsWith(filename, ".csv")) {
      safe_write_file(
        data = metadata,
        path = output_file_path,
        write_fn = write.csv,
        verbose = RUNTIME_CONFIG$debug_verbose,
        interactive = interactive(),  # Use interactive() to detect if running interactively
        row.names = FALSE
      )
    } else if (endsWith(filename, ".tsv")) {
      safe_write_file(
        data = bmc_metadata,
        path = output_file_path,
        write_fn = write.table,
        verbose = RUNTIME_CONFIG$debug_verbose,
        interactive = interactive(),
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    } else if (endsWith(filename, ".R")) {
      safe_write_file(
        data = EXPERIMENT_CONFIGURATION_PATH,
        path = output_file_path,
        write_fn = file.copy,
        verbose = RUNTIME_CONFIG$debug_verbose,
        interactive = interactive(),
        overwrite = TRUE
      )
    } else {
      warning(sprintf("Unsupported file extension for file: %s", filename))
    }
  }
}
cat("Final metadata organization\n")
cat("----------------------\n")
print(metadata)
cat("----------------------\n")
