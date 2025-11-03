#!/usr/bin/env Rscript
################################################################################
# BMC ChIP-seq Experiment Setup
# Author: Luis | Date: 2025-11-02 | Version: 3.0.0
################################################################################
# PURPOSE:
#   Creates standardized directory structure and metadata files for BMC ChIP-seq experiments
#
# USAGE:
#   1) Update configuration_experiment_bmc.R file.
#   2) From R REPL: source("core_scripts/setup_bmc_experiment.R")
#
# DEPENDENCIES: 
#   ./core_scripts/configuration_experiment_bmc.R
#   ./core_scripts/configuration_script_bmc.R
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

category_names <- names(EXPERIMENT_CONFIG$CATEGORIES)
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

# Validate categories - simple for loop
for (category_name in names(EXPERIMENT_CONFIG$CATEGORIES)) {
  category_values <- EXPERIMENT_CONFIG$CATEGORIES[[category_name]]

  if (!is.character(category_values)) {
    stop(sprintf("Category '%s' must be character vector", category_name))
  }

  if (any(duplicated(category_values))) {
    duplicated_values <- category_values[duplicated(category_values)]
    stop(sprintf("Category '%s' has duplicates: %s",
                category_name, paste(duplicated_values, collapse = ", ")))
  }
}

# Validate invalid combinations - simple for loop
for (combo_name in names(EXPERIMENT_CONFIG$INVALID_COMBINATIONS)) {
  combo_expr <- EXPERIMENT_CONFIG$INVALID_COMBINATIONS[[combo_name]]
  referenced_vars <- all.vars(combo_expr)
  invalid_refs <- setdiff(referenced_vars, category_names)

  if (length(invalid_refs) > 0) {
    stop(sprintf("INVALID_COMBINATIONS '%s' references unknown categories: %s",
                combo_name, paste(invalid_refs, collapse = ", ")))
  }
}


# Validate COLUMN_ORDER matches CATEGORIES
if (!identical(sort(names(EXPERIMENT_CONFIG$CATEGORIES)),
               sort(EXPERIMENT_CONFIG$COLUMN_ORDER))) {
  stop("COLUMN_ORDER must include exactly all category names")
}

# Validate CONTROL_FACTORS references valid categories
for (factor_name in names(EXPERIMENT_CONFIG$CONTROL_FACTORS)) {
  factor_categories <- EXPERIMENT_CONFIG$CONTROL_FACTORS[[factor_name]]
  invalid_refs <- setdiff(factor_categories, names(EXPERIMENT_CONFIG$CATEGORIES))

  if (length(invalid_refs) > 0) {
    stop(sprintf("CONTROL_FACTORS '%s' references unknown categories: %s",
                factor_name, paste(invalid_refs, collapse = ", ")))
  }
}

#====================
# Source script configuration file
#====================
source(SCRIPT_CONFIGURATION_PATH)

stopifnot(
  "EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID does not match configuration EXPERIMENT_ID" =
    identical(EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID, EXPERIMENT_IDS),
  "Experiment configuration file does not exist. Please copy the template." =
    file.exists(EXPERIMENT_CONFIGURATION_PATH),
  "Experiment ID must be a character string" =
    is.character(EXPERIMENT_IDS),
  "Invalid experiment ID format. Expected: YYMMDD'Bel'" =
    grepl("^\\d{6}Bel$", EXPERIMENT_IDS),
  "Script EXPERIMENT_IDS is not the same as CONFIG EXPERIMENT_IDS" =
    EXPERIMENT_IDS == EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID,
  "Only one experiment id required for this script" =
    length(EXPERIMENT_IDS) == 1
)

message("Configuration file sourced... Checking configuration variables.")

# Ensure the variables expected in the script were //
# defined in the configuration file. //
# configuration_script_bmc.R //
required_configuration_variables <- c(
  "EXPERIMENT_IDS", "EXPERIMENT_DIR",
  "DATA_DIRECTORIES"
)

is_missing_variable <- !sapply(required_configuration_variables, exists)
missing_variables <- required_configuration_variables[is_missing_variable]

if (length(missing_variables) > 0 ) {
  stop("Missing variable. Please define in 'configuration_script_bmc.R' file.",
     paste(missing_variables, collapse = ", "))
}
message("All variables defined in the configuration file...")

#-------------------------------------------------------------------------------
# Directory Setup
#-------------------------------------------------------------------------------
# Create directory structure
full_paths <- file.path(EXPERIMENT_DIR, DATA_DIRECTORIES)
for (path in full_paths) {
  if (RUNTIME_CONFIG$output_dry_run) {
    cat(sprintf("[DRY RUN] Would create: %s\n", path))
    next

  }

  dir_created <- dir.create(path, recursive = TRUE, showWarnings = FALSE)

  if (RUNTIME_CONFIG$debug_verbose) {
    cat(sprintf(
      "[%s] %s\n", 
      if (dir_created) "CREATED" else "EXISTS", 
      path
    ))

  }

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
# Use do.call for multi-column ordering
metadata <- metadata[
  do.call(
    order,
    metadata[EXPERIMENT_CONFIG$COLUMN_ORDER]
  ),
]

# Verify sample count
n_samples <- nrow(metadata)
expected <- EXPERIMENT_CONFIG$METADATA$EXPECTED_SAMPLES_SETUP
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
# Apply row order correction if specified in config
if (!is.null(EXPERIMENT_CONFIG$ROW_ORDER_CORRECTION)) {
  from_row <- EXPERIMENT_CONFIG$ROW_ORDER_CORRECTION$from_row
  to_row <- EXPERIMENT_CONFIG$ROW_ORDER_CORRECTION$to_row

  # Validate inputs
  stopifnot(
    is.numeric(from_row),
    length(from_row) == 1,
    is.numeric(to_row),
    length(to_row) == 1,
    from_row >= 1,
    from_row <= nrow(metadata),
    to_row >= 1,
    to_row <= nrow(metadata),
    from_row != to_row
  )

  # Start with the original row indices (1, 2, 3, ..., N)
  original_row_indices <- seq_len(nrow(metadata))

  # Remove the row that will be moved, creating a shortened list of row indices
  row_indices_without_moved_sample <- original_row_indices[original_row_indices != from_row]

  # Determine the insertion point in the shortened list so that the moved sample
  # ends up at position 'to_row' in the final reordered metadata
  insertion_index_in_shortened_list <- to_row - 1

  # Insert the moved row at the calculated position to produce the final row order
  final_row_order <- append(
    x     = row_indices_without_moved_sample,
    values = from_row,
    after = insertion_index_in_shortened_list
  )

  metadata <- metadata[final_row_order, , drop = FALSE]

  moved_sample_final_position <- which(final_row_order == from_row)
  cat("Moved sample (original row", from_row, ") is now at position:", moved_sample_final_position, "\n")
  stopifnot(moved_sample_final_position == to_row)
}

#-------------------------------------------------------------------------------
# Remove samples that dropped out
#-------------------------------------------------------------------------------
# Remove sample that was dropped during preparation (if specified)
if (!is.null(EXPERIMENT_CONFIG$SAMPLE_DROPOUT_CONDITION)) {
  # Evaluate the quoted expression in the context of the metadata dataframe
  dropout_rows <- eval(EXPERIMENT_CONFIG$SAMPLE_DROPOUT_CONDITION, envir = metadata)
  #dropout_rows <- length(EXPERIMENT_CONFIG$SAMPLE_DROPOUT_CONDITION)

  # Validate: must be logical and same length as rows
  stopifnot(
    is.logical(dropout_rows),
    length(dropout_rows) == nrow(metadata)
  )

  # Ensure exactly one sample is being dropped (adjust if you ever drop multiple)
  # Replace with length or other way to measure amount to drop
  stopifnot(sum(dropout_rows, na.rm = TRUE) == 1)

  # Remove the dropped sample
  metadata <- metadata[!dropout_rows, , drop = FALSE]
}

# Verify sample count after processing
n_samples <- nrow(metadata)
expected <- EXPERIMENT_CONFIG$METADATA$EXPECTED_SAMPLES_POST
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
# Define output files
output_files <- list(
  sample_grid = list(
    data = metadata_grid,
    path = file.path(EXPERIMENT_DIR, "documentation", sprintf("%s_sample_grid.csv", EXPERIMENT_IDS)),
    type = "csv"
  ),
  bmc_table = list(
    data = bmc_table,
    path = file.path(EXPERIMENT_DIR, "documentation", sprintf("%s_bmc_table.tsv", EXPERIMENT_IDS)),
    type = "tsv"
  ),
  experiment_configuration = list(
    data = file.path(CORE_SCRIPTS_PATH, "configuration_experiment_bmc.R"),
    path = file.path(EXPERIMENT_DIR, "documentation", sprintf("%s_configuration_experiment_bmc.R", EXPERIMENT_IDS)),
    type = "r"
  )
)

# Write files
for (file_name in names(output_files)) {
  # Extract to descriptive variables
  file_info <- output_files[[file_name]]
  output_data <- file_info$data
  output_path <- file_info$path
  file_type <- file_info$type

  if (OUTPUT_DRY_RUN) {
    status <- if (file.exists(output_path)) "OVERWRITE" else "CREATE"
    message(sprintf("[DRY RUN] Would %s: %s", status, output_path))
    next
  }

  # Check for existing file in interactive mode
  if (interactive() && file.exists(output_path)) {
    user_input <- readline(prompt = sprintf("File exists: %s. Overwrite? (y/n): ",
                                            output_path))
    if (tolower(user_input) != "y") {
      message(sprintf("Skipped: %s", output_path))
      next
    }
  }

  # Write based on file type
  tryCatch({
    switch(file_type,
      csv = write.csv(output_data, output_path, row.names = FALSE),
      tsv = write.table(output_data, output_path, sep = "\t",
                        row.names = FALSE, quote = FALSE),
      r = file.copy(output_data, output_path, overwrite = TRUE),
      rds = saveRDS(output_data, output_path),
      stop(sprintf("Unknown file type: %s", file_type))
    )
    message(sprintf("Wrote: %s", output_path))

  }, error = function(e) {
    warning(sprintf("Failed to write %s: %s", output_path, e$message))

  })

}
