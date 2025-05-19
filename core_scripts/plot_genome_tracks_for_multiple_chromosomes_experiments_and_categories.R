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
missing_variables <- required_configuration_variables[!sapply(required_configuration_variables, exists)]
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
# Refactor
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
metadata_aligned <- lapply(metadata_list,
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
# Oh. Can this be removed to the configuration file as well?
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
  stop(sprintf(
    fmt = "No reference genome files found matching pattern '%s' in: %s",
    FILE_GENOME_DIRECTORY,
    FILE_FEATURE_DIRECTORY
  ))
}

if (!file.exists(REF_GENOME_FILE)) {
  stop(sprintf("Reference genome file not accessible: %s", REF_GENOME_FILE[1]))
}
if (length(FEATURE_FILE) == 0) {
  warning(sprintf(
    fmt = "No feature files found matching pattern '%s' in: %s",
    FILE_FEATURE_PATTERN,
    FILE_FEATURE_DIRECTORY
  ))
}

REFERENCE_GENOME_DSS <- Biostrings::readDNAStringSet(REF_GENOME_FILE)
CHROMOSOME_WIDTHS <- REFERENCE_GENOME_DSS[CHROMOSOMES_TO_PLOT]@ranges@width
CHROMOSOMES_IN_ROMAN <- paste0("chr", utils::as.roman(CHROMOSOMES_TO_PLOT))

GENOME_RANGE_TO_LOAD <- GenomicRanges::GRanges(
    seqnames = CHROMOSOMES_IN_ROMAN,
    ranges = IRanges::IRanges(start = 1, end = CHROMOSOME_WIDTHS),
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
  # GENOME_FEATURES <- GENOME_FEATURES[seqnames(GENOME_FEATURES) == CHROMOSOMES_IN_ROMAN]
  GENOME_FEATURES <- GenomeInfoDb::keepSeqlevels(GENOME_FEATURES, CHROMOSOMES_IN_ROMAN, pruning.mode = "coarse")
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
# TODO: Move to configuration. //
row_filtering_expression <- quote(!(rescue_allele == "4R" & suppressor_allele == "NONE"))
if (exists("row_filtering_expression") & !is.null(row_filtering_expression)) {
  expr_vars <- all.vars(row_filtering_expression)
  missing_expr_vars <- setdiff(expr_vars, colnames(metadata_df))
  if (length(missing_expr_vars) > 0) {
    stop("Expression refers to columns that don't exist in metadata_df: ",
         paste(missing_expr_vars, collapse = ", "))
  }
  tryCatch({
    # would need to add more conditionals if I assign to a new variable //
    metadata_df <- metadata_df[eval(row_filtering_expression, envir = metadata_df), ]
  }, error = function(e) {
    stop("Error evaluating row_filtering_expression: ", e$message)
  })
} # end if error handling for expression

# Define columns to compare by and exclude columns for grouping
# TODO: Move to configuration. //
target_comparison_columns <- c("rescue_allele", "suppressor_allele")
plot_name_comparison_column_section <- paste(
  gsub("_", "", target_comparison_columns),
  collapse = "."
)
metadata_columns_to_exclude <- c(
    "sample_type", "sample_ids",
    "bigwig_file_paths", "full_name",
    "short_name"
)

missing_excluded_columns <- setdiff(metadata_columns_to_exclude, colnames(metadata_df))
missing_comparison_columns <- setdiff(target_comparison_columns, colnames(metadata_df))

# Validate target comparison columns exist
if (length(missing_comparison_columns) > 0) {
  stop("Missing required comparison columns in metadata_df: ",
       paste(missing_comparison_columns, collapse = ", "))
}
# Validate excluded columns (warning only, not critical)
if (length(missing_excluded_columns) > 0) {
  warning("Some excluded columns don't exist in metadata_df: ",
          paste(missing_excluded_columns, collapse = ", "))
  # Remove non-existent columns from the exclusion list
  metadata_columns_to_exclude <- intersect(metadata_columns_to_exclude, colnames(metadata_df))
}

experimental_condition_columns <- setdiff(
  colnames(metadata_df),
  union(target_comparison_columns, metadata_columns_to_exclude)
)
metadata_df$track_name <- do.call(paste,
  args = c(
    lapply(metadata_df[, target_comparison_columns],
           as.character),
    sep = "-"
  )
)
experimental_condition_data <- metadata_df[, experimental_condition_columns, drop = FALSE]
metadata_df$experimental_condition_id <- do.call(
  paste,
  args = c(experimental_condition_data, sep = "|")
)

unique_experimental_conditions <- unique(metadata_df$experimental_condition_id)
unique_condition_data <- unique(experimental_condition_data)

unique_condition_text_df <- data.frame(
  lapply(unique_condition_data, as.character),
  stringsAsFactors = FALSE
)
# For each row, create "column_name: value" strings
experimental_condition_titles <- apply(
  unique_condition_text_df, 1,
  function(row_values) {
    paste(names(row_values),
      row_values,
      sep = ": ",
      collapse = "\n"
    )
  }
)

total_number_of_conditions <- length(unique_experimental_conditions)
total_number_of_samples <- nrow(metadata_df)
total_number_of_chromosomes <- length(CHROMOSOMES_TO_PLOT)

stopifnot(
  "No experimental condition columns remain. See experimental_condition_columns" =
    length(experimental_condition_columns) > 0,
  "No unique experimental conditions identified. See unique_experimental_conditions." =
    total_number_of_conditions > 0,
  "No samples in metadata df. See metadata_df." =
    total_number_of_samples > 0,
  "No chromosomes to plot. See CHROMOSOMES_TO_PLOT." =
    total_number_of_chromosomes > 0,
  "Expect same number of titles and conditions. See condition_plot_titles and unique_experimental_condtions." =
    length(experimental_condition_titles) == total_number_of_conditions
)

for (condition_idx in seq_len(total_number_of_conditions)) {
  message("=== For loop for group ===")
  message(sprintf(
    fmt = "  Processing group: %s / %s ",
    condition_idx, total_number_of_conditions)
  )
  # Iteration setup ----------
  current_condition <- unique_experimental_conditions[condition_idx]
  current_condition_base_title <- experimental_condition_titles[condition_idx]
  is_condition_row <- metadata_df$experimental_condition_id == current_condition
  current_condition_df <- metadata_df[is_condition_row, ]
  current_number_of_samples <- nrow(current_condition_df)
  # Move to each or inside the most nested
  debug_print(list(
    "title" = "Debug group plotting",
    ".Number of rows" = current_number_of_samples,
    ".Number of original rows" = total_number_of_samples,
    ".Current condition" = current_condition,
    ".Current title" = current_condition_base_title
  ))

  if (current_number_of_samples == 0) {
    warning(sprintf(
      fmt = "Condition '%s' does not have any samples. ",
      current_condition
    ))
    next
  }

  for (chromosome_idx in seq_len(total_number_of_chromosomes)) {
    message("  --- For loop for chromosome ---")
    message(sprintf(
      fmt = "    Processing chromosome: %s / %s ",
      chromosome_idx, total_number_of_chromosomes
    ))
    current_chromosome <- CHROMOSOMES_IN_ROMAN[chromosome_idx]
    chromosome_title_section <- paste0("Chromosome: ", current_chromosome)
    #current_genome_range_to_load <- GENOME_RANGE_TO_LOAD[chromosome_idx]
    current_condition_title <- paste(
      current_condition_base_title,
      chromosome_title_section,
      sep = "\n"
    )
    track_container <- list(
      Gviz::GenomeAxisTrack(
        name = sprintf("Chr %s Axis", current_chromosome),
        fontcolor.title = "black",
        cex.title = 0.7,
        background.title = "white"
      )
    )
    plot_file_name <- paste(
      "condition_plots_",
      current_chromosome,
      gsub("|", ".", current_condition),
      "_",
      plot_name_comparison_column_section,
      ".svg",
      sep = ""
    )
  debug_print(list(
    "title" = "Debug chromosome",
    ".Current chromosome" = current_chromosome,
    ".Plot name" = plot_file_name
  ))
    # Move the sample for loop here.
    message("  --- end chromosome iteration ---")
  }
    for (sample_idx in seq_len(current_number_of_samples)) {
      message("    --- For loop for sample ---")
      message(sprintf(
        fmt = "      Processing row: %s / %s ",
        sample_idx, current_number_of_samples
      ))
      current_sample_track_name <- current_condition_df$track_name[sample_idx]
      current_bigwig_file_path <- current_condition_df$bigwig_file_paths[sample_idx]
      container_length <- length(track_container)
      debug_print(list(
        "title" = "Sample iteration",
        ".Track_name" = current_sample_track_name,
        ".Bigwig file path" = current_bigwig_file_path,
        ".Container Length" = container_length
      ))
      #bigwig_data <- rtracklayer::import(
      #  current_bigwig_file_path,
      #  format = "BigWig",
      #  which = current_genome_range_to_load
      #)
      #track_container[[container_length + 1]] <- Gviz::DataTrack(
      #    range = bigwig_data,
      #    name = track_name,
      #    # Apply styling
      #    showAxis = TRUE,
      #    showTitle = TRUE,
      #    type = "h",
      #    size = 1.2,
      #    background.title = "white",
      #    fontcolor.title = "black",
      #    col.border.title = "#e0e0e0",
      #    cex.title = 0.7,
      #    fontface = 1,
      #    title.width = 1.0,
      #    col = "darkblue",
      #    fill = "darkblue"
      #)
      # Grab the appropriate data. Load the data
      message("    --- end row iteration ---")
    }

  # [4] Annotation Track (Conditional)
  #if (exists("GENOME_FEATURES")) {
  #  track_container[[length(track_container)]] <- Gviz::AnnotationTrack(
  #    GENOME_FEATURES,
  #    name = "Features",
  #    size = 0.5,
  #    background.title = "lightgray",
  #    fontcolor.title = "black",
  #    showAxis = FALSE,
  #    background.panel = "#f5f5f5",
  #    cex.title = 0.6,
  #    fill = "#8b4513",
  #    col = "#8b4513"
  #  )
  #}

  message("  Ended condition iteration ", condition_idx)
  message("=== end condition iteration ===")
  message("\n")
}
