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
# See template_script_configuration.R or script_configuration.R
required_configuration_variables <- c(
  "EXPERIMENT_IDS",
  "EXPERIMENT_DIR",
  "CHROMOSOMES_TO_PLOT",
  "OUTPUT_FORMAT",
  "ACCEPT_CONFIGURATION",
  "SKIP_PACKAGE_CHECKS"
)
missing_variables <- required_configuration_variables[!exists(required_configuration_variables)]
if (length(missing_variables) > 0 ) {
  stop("Missing variable. Please define in 'script_configuration.R' file.",
       paste(missing_variables, collapse = ", "))
}
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
#bigwig_pattern <- "processed_.*_sequence_to_S288C_blFiltered_CPM\\.bw$"
bigwig_pattern <- "processed_.*_sequence_to_S288C_blFiltered_CPM\\.bw$"
fastq_pattern <- "consolidated_.*_sequence\\.fastq$"
sample_id_capture_pattern <- "consolidated_([0-9]{1,6})_sequence\\.fastq$"

expected_number_of_samples <- 0
REQUIRED_DIRECTORIES <- c("fastq", "coverage")
for (experiment_idx in seq_len(number_of_experiments)) {
  message("--- Start experiment idx ---")
  current_experiment_path <- EXPERIMENT_DIR[experiment_idx]
  current_config_path <- config_paths[experiment_idx]
  current_metadata_path <- metadata_paths[experiment_idx]
  # Build all required directory paths for this experiment
  required_data_paths <- file.path(current_experiment_path, REQUIRED_DIRECTORIES)
  names(required_data_paths) <- REQUIRED_DIRECTORIES
  missing_dirs <- required_data_paths[!dir.exists(required_data_paths)]
  if (length(missing_dirs) > 0) {
    stop("Missing required experiment subdirectories: ",
         paste(missing_dirs, collapse = ", "))
  }
  debug_print(list(
    "title" = "Debug Path information",
    ".Current Experiment path" = current_experiment_path,
    ".Current metadata path" = current_metadata_path,
    ".Current Config path" = current_config_path,
    ".Fastq path" = required_data_paths[["fastq"]],
    ".Coverage path" = required_data_paths[["coverage"]]
  ))

  # Load *_CONFIG variables ---------
  source(current_config_path)
  current_metadata_df <- read.csv(current_metadata_path, stringsAsFactors = FALSE)

  # Gather additional metadata --------
  current_metadata_df$experiment_id <- EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID
  # Find fastq files and extract sample IDs
  fastq_files <- list.files(
    path = required_data_paths[["fastq"]],
    pattern = fastq_pattern,
    full.names = FALSE
  )
  bigwig_files <- list.files(
    path = required_data_paths[["coverage"]],
    pattern = bigwig_pattern,
    full.names = TRUE
  )
  #bigwig_basenames <- basename(bigwig_files)
  stopifnot(
    "No fastq files found." = length(fastq_files) > 0,
    "No bigwig files found." = length(bigwig_files) > 0
    )
  sample_ids <- gsub(
    pattern = sample_id_capture_pattern,
    replacement = "\\1",
    x = fastq_files
  )

  stopifnot(
     "Length of samples_ids is not lower than length of fastq files." =
        all(nchar(sample_ids) < nchar(fastq_files)),
     "Some of the extracted sample_ids are empty." =
        all(nchar(sample_ids) > 0),
     "Extracted sample ids cannot be coerced to numeric." =
       all(!is.na(as.numeric(sample_ids))),
     "Number of extracted sample_ids does not match number of current metadata df rows." =
       nrow(current_metadata_df) == length(sample_ids)
  )

  expected_number_of_samples <- expected_number_of_samples + EXPERIMENT_CONFIG$METADATA$EXPECTED_SAMPLES
  metadata_categories_list[[experiment_idx]] <- EXPERIMENT_CONFIG$CATEGORIES
  current_metadata_df$bigwig_file_paths <- bigwig_files
  current_metadata_df$sample_ids <- sample_ids
  debug_print(list(
    "title" = "Debug metadata loading",
    ".Number of rows" = nrow(current_metadata_df),
    ".Number of samples ids" = length(sample_ids),
    ".Number of bigwig files" = length(bigwig_files),
    ".Number of expected samples" = EXPERIMENT_CONFIG$METADATA$EXPECTED_SAMPLES
  ))
  # Add determination of sample ids and addition to metadata frame
  metadata_list[[experiment_idx]] <- current_metadata_df
  message("--- End iteration ---")
  message("\n")
} # end for loop
message("Finished metadata loading...")

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
stopifnot("Metadata does not have expected number of rows." = nrow(metadata_df) == expected_number_of_samples)

# Convert the columns of the metadata_df to factors.
for (col_name in intersect(names(merged_categories), colnames(metadata_df))) {
  metadata_df[[col_name]] <- factor(
    metadata_df[[col_name]],
    levels = merged_categories[[col_name]],
    ordered = TRUE
    )
}
message("Finished metadata processing...")
#metadata_df <- metadata_df[do.call(order, metadata_df[intersect(EXPERIMENT_CONFIG$COLUMN_ORDER, colnames(metadata_df))]), ]
################################################################################
# Setup genome and feature files
################################################################################
# Ensure supplementary files for plotting genome tracks are present
FILE_GENOME_DIRECTORY <- file.path(Sys.getenv("HOME"), "data", "REFGENS")
FILE_GENOME_PATTERN <- "S288C_refgenome.fna"
FILE_FEATURE_DIRECTORY <- file.path(Sys.getenv("HOME"), "data", "feature_files")
FILE_FEATURE_PATTERN <- "eaton_peaks"
stopifnot(
    "Genome directory not found" = dir.exists(FILE_GENOME_DIRECTORY),
    "Feature directory not found" = dir.exists(FILE_FEATURE_DIRECTORY)
)

# Load reference genome
REF_GENOME_FILE <- list.files(
    path = FILE_GENOME_DIRECTORY,
    pattern = FILE_GENOME_PATTERN,
    full.names = TRUE,
    recursive = TRUE
)[1]

# Load feature file (annotation)
FEATURE_FILE <- list.files(
    path = FILE_FEATURE_DIRECTORY,
    pattern = FILE_FEATURE_PATTERN,
    full.names = TRUE
)[1]

if (length(REF_GENOME_FILE) == 0) {
    stop(sprintf("No reference genome files found matching pattern '%s' in: %s",
                FILE_GENOME_DIRECTORY,
                FILE_FEATURE_DIRECTORY))
}

if (!file.exists(REF_GENOME_FILE)) {
    stop(sprintf("Reference genome file not accessible: %s", REF_GENOME_FILE[1]))
}
if (length(FEATURE_FILE) == 0) {
    warning(sprintf("No feature files found matching pattern '%s' in: %s",
                   FILE_FEATURE_PATTERN,
                   FILE_FEATURE_DIRECTORY))
}

REFERENCE_GENOME_DSS <- Biostrings::readDNAStringSet(REF_GENOME_FILE)
CHROMOSOME_WIDTH <- REFERENCE_GENOME_DSS[CHROMOSOME_TO_PLOT]@ranges@width
CHROMOSOME_ROMAN <- paste0("chr", utils::as.roman(CHROMOSOME_TO_PLOT))

GENOME_RANGE_TO_LOAD <- GenomicRanges::GRanges(
    seqnames = CHROMOSOME_ROMAN,
    ranges = IRanges::IRanges(start = 1, end = CHROMOSOME_WIDTH),
    strand = "*"
)

if (!is.null(FEATURE_FILE)) {
  GENOME_FEATURES <- rtracklayer::import(FEATURE_FILE)
  # Convert to chrRoman format
  GenomeInfoDb::seqlevels(GENOME_FEATURES) <- paste0(
    "chr",
    utils::as.roman(gsub("chr", "", GenomeInfoDb::seqlevels(GENOME_FEATURES)))
  )
  # Fitler out the rest of the chromosomes
  # GENOME_FEATURES <- GENOME_FEATURES[seqnames(GENOME_FEATURES) == CHROMOSOME_ROMAN]
  GENOME_FEATURES <- GenomeInfoDb::keepSeqlevels(GENOME_FEATURES, CHROMOSOME_ROMAN, pruning.mode = "coarse")
}

####################
# Plot bigwig files
# For each repeat,
#    subset by control column,
#    plot the tracks for multiple chromosomes
####################
# Define any rows to exclude
# For example, remove rows that you dont want to be included in the plots
# Prefilter the metadata_df
row_filtering_expression <- quote(rescue_allele == "4R" & suppressor_allele == "NONE")

target_comparison_columns <- c("rescue_allele", "suppressor_allele")
metadata_columns_to_exclude <- c(
    "sample_type", "sample_ids",
    "bigwig_file_paths", "full_name",
    "short_name"
)
experimental_condition_columns <- setdiff(
  colnames(metadata_df),
  union(target_comparison_columns, metadata_columns_to_exclude)
)
metadata_df$experimental_condition_id <- do.call(
  paste,
  c(metadata_df[experiment_condition_columns], sep = "|")
)
unique_experimental_conditions <- unique(metadata_df$experimental_condition_id)
total_number_of_conditions <- length(unique_experimental_conditions)


for (condition_idx in seq_len(total_number_of_conditions)) {
  message("=== For loop for group ===")
  message(sprintf(
    fmt = "  Processing group: %s / %s ",
    condition_idx, total_number_of_conditions)
  )
  current_condition <- unique_experimental_conditions[condition_idx]
  is_condition_row <- metadata_df$experiment_condition_id == current_condition
  current_condition_df <- metadata_df[is_condition_row, ]
  number_of_samples <- nrow(current_condition_df)
  debug_print(list(
    "title" = "Debug group plotting",
    ".Number of rows" = number_of_samples,
    ".Current group" = current_condition
  ))
  if (number_of_samples == 0) {
    warning(sprintf(
      fmt = "Condition '%s' does not have any samples. ",
      current_condition
    ))
    next
  }

  for (sample_idx in seq_len(number_of_samples)) {
    message("  --- For loop for row ---")
    message(sprintf(
      fmt = "    Processing row: %s / %s ",
      sample_idx, number_of_samples
    ))
    # Grab the appropriate data. Load the data
    message("  --- end row iteration ---")
  }

  message(sprintf(
    fmt = "  Ended iteration %s",
    condition_idx
  ))
  message("=== end group iteration ===")
  message("\n")
}
#DF_EXPERIMENT_IDS <- unique(metadata_df$experiment_id)
#NUMBER_OF_EXPERIMENT_IDS <- length(DF_EXPERIMENT_IDS)
#for (experiment_idx in seq_len(NUMBER_OF_EXPERIMENT_IDS)) {
  #message("--- For loop for metadata ---")
  #message(sprintf(
    #fmt = "  Processing row: %s / %s ",
    #experiment_idx, NUMBER_OF_EXPERIMENT_IDS)
  #)
  #current_experiment_id <- DF_EXPERIMENT_IDS[experiment_idx]
  #is_experiment_id_row <- metadata_df$experiment_id == current_experiment_id
  #experiment_id_subset_df <- metadata_df[is_experiment_id_row, ]
#}
