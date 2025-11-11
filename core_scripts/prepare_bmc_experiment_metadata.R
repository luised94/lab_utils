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

required_packages <- c(
  "rtracklayer", "GenomicRanges",
  "Biostrings", "RColorBrewer"
)

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

#===============================================================================
# SECTION 5: Discover and Match Files
#===============================================================================
if (RECOMPUTATION_FLAGS$FILE_PATHS) {
  message("=== SECTION 5: Discover and Match Files ===")

  # Define data directories
  fastq_dir <- file.path(EXPERIMENT_DIR_PATH, "fastq")
  coverage_dir <- file.path(EXPERIMENT_DIR_PATH, "coverage")

  # Validate directories exist
  if (!dir.exists(fastq_dir)) {
    stop("FASTQ directory not found: ", fastq_dir, "\n",
         "Expected: <experiment>/fastq/processed/")
  }

  if (!dir.exists(coverage_dir)) {
    stop("Coverage directory not found: ", coverage_dir, "\n",
         "Expected: <experiment>/coverage/")
  }

  # Find FASTQ files
  fastq_files <- list.files(
    path = fastq_dir,
    pattern = FASTQ_PATTERN,
    full.names = FALSE
  )

  if (length(fastq_files) == 0) {
    stop("No FASTQ files found in: ", fastq_dir, "\n",
         "Pattern: ", FASTQ_PATTERN)
  }

  message("Found ", length(fastq_files), " FASTQ files")

  # Extract sample IDs from FASTQ filenames
  sample_ids <- gsub(SAMPLE_ID_PATTERN, "\\1", fastq_files)
  sample_ids <- unique(sub(".*-", "", sample_ids))

  # Validate extraction succeeded
  if (any(nchar(sample_ids) == 0)) {
    stop("Failed to extract sample IDs from FASTQ files.\n",
         "Check SAMPLE_ID_PATTERN: ", SAMPLE_ID_PATTERN)
  }

  # Validate sample ID count matches metadata
  if (length(sample_ids) != nrow(metadata_df)) {
    stop("Sample ID count mismatch.\n",
         "FASTQ files: ", length(sample_ids), "\n",
         "Metadata rows: ", nrow(metadata_df), "\n",
         "Sample IDs: ", paste(sample_ids, collapse = ", "))
  }

  # Add sample IDs to metadata
  metadata_df$sample_ids <- sample_ids
  message("Sample IDs added to metadata")

  # Find bigwig files
  bigwig_files <- list.files(
    path = coverage_dir,
    pattern = BIGWIG_PATTERN,
    full.names = TRUE
  )

  if (length(bigwig_files) == 0) {
    warning("No bigwig files found in: ", coverage_dir, "\n",
            "Pattern: ", BIGWIG_PATTERN, "\n",
            "Bigwig paths will be NA")
    metadata_df$bigwig_file_path <- NA

  } else {
    message("Found ", length(bigwig_files), " bigwig files")

    # Extract sample IDs from bigwig filenames
    bigwig_basenames <- basename(bigwig_files)
    bigwig_sample_ids <- sub(".*_(D[0-9]{2}-[0-9]{1,6})_.*", "\\1", bigwig_basenames)

    # Create lookup table
    # Works even if some bigwig files are missing.
    bigwig_lookup <- setNames(bigwig_files, bigwig_sample_ids)

    # Match bigwig files to metadata by sample ID
    metadata_df$bigwig_file_path <- bigwig_lookup[metadata_df$sample_ids]

    # Report matching statistics
    missing_bigwig <- is.na(metadata_df$bigwig_file_path)
    n_matched <- sum(!missing_bigwig)
    n_missing <- sum(missing_bigwig)

    message("Matched bigwig files: ", n_matched, " of ", nrow(metadata_df))

    if (n_missing > 0) {
      missing_samples <- metadata_df$sample_ids[missing_bigwig]
      warning("Samples missing bigwig files:\n",
              paste(missing_samples, collapse = ", "))
    }
  }

  message("Section 5 complete\n")

} else {
  message("=== SECTION 5: Skipped (FILE_PATHS = FALSE) ===\n")

}


#===============================================================================
# SECTION 6: Compute Track Limits
#===============================================================================
if (RECOMPUTATION_FLAGS$TRACK_LIMITS) {
  message("=== SECTION 6: Compute Track Limits ===")
  
  # Check if bigwig paths exist in metadata
  if (!"bigwig_file_path" %in% colnames(metadata_df)) {
    stop("Column 'bigwig_file_path' not found in metadata.\n",
         "Run Section 5 (FILE_PATHS) first")
  }
  
  # Initialize track limit columns
  metadata_df$track_ymin <- 0
  metadata_df$track_ymax <- 100  # Default
  
  # Define chromosomes to analyze
  chromosomes_roman <- paste0("chr", utils::as.roman(CHROMOSOMES_TO_PLOT))
  
  # Count samples with valid bigwig files
  has_bigwig <- !is.na(metadata_df$bigwig_file_path) & 
                file.exists(metadata_df$bigwig_file_path)
  n_with_bigwig <- sum(has_bigwig)
  
  if (n_with_bigwig == 0) {
    warning("No valid bigwig files to compute track limits.\n",
            "Using default limits: 0-100")
    message("Section 6 complete\n")
  } else {
    message("Computing track limits for ", n_with_bigwig, " samples...")
    message("Chromosomes: ", paste(CHROMOSOMES_TO_PLOT, collapse = ", "))
    
    n_processed <- 0
    
    for (row_idx in seq_len(nrow(metadata_df))) {
      bigwig_path <- metadata_df$bigwig_file_path[row_idx]
      
      # Skip if no bigwig file
      if (is.na(bigwig_path) || !file.exists(bigwig_path)) {
        next
      }
      
      # Import and process bigwig data
      tryCatch({
        # Load entire bigwig
        bigwig_gr <- rtracklayer::import(bigwig_path)
        
        # Filter to chromosomes of interest
        bigwig_gr <- bigwig_gr[seqnames(bigwig_gr) %in% chromosomes_roman]
        
        if (length(bigwig_gr) > 0) {
          # Get score values
          scores <- bigwig_gr$score
          scores <- scores[!is.na(scores) & is.finite(scores)]
          
          if (length(scores) > 0) {
            metadata_df$track_ymin[row_idx] <- 0
            metadata_df$track_ymax[row_idx] <- max(scores, na.rm = TRUE)
          }
        }
        
        n_processed <- n_processed + 1
        
        # Progress message every 10 samples
        if (n_processed %% 10 == 0) {
          message("  Processed ", n_processed, " of ", n_with_bigwig)
        }
        
      }, error = function(e) {
        warning("Failed to process bigwig: ", basename(bigwig_path), "\n",
                "Error: ", e$message)
      })
    }
    
    message("Track limits computed for ", n_processed, " samples")
    message("Section 6 complete\n")
  }
  
} else {
  message("=== SECTION 6: Skipped (TRACK_LIMITS = FALSE) ===\n")
}

#===============================================================================
# SECTION 7: Generate Color Scheme
#===============================================================================
if (RECOMPUTATION_FLAGS$COLOR_SCHEMES) {
  message("=== SECTION 7: Generate Color Scheme ===")

  # Define category to color by (make this configurable later)
  category_to_color_by <- "antibody"

  # Validate column exists
  if (!category_to_color_by %in% colnames(metadata_df)) {
    stop("Category not found in metadata: ", category_to_color_by, "\n",
         "Available columns: ", paste(colnames(metadata_df), collapse = ", "))
  }

  # Get unique category values
  category_values <- unique(metadata_df[[category_to_color_by]])
  number_of_colors <- length(category_values)

  message("Generating colors for: ", category_to_color_by)
  message("  Unique values: ", number_of_colors)

  # Generate colors with reproducible seed
  category_seed <- sum(utf8ToInt(category_to_color_by))
  set.seed(category_seed)

  if (number_of_colors <= 8) {
    # Use ColorBrewer palette for small sets
    category_palette <- RColorBrewer::brewer.pal(
      max(3, number_of_colors),
      "Set2"
    )
    category_colors <- category_palette[seq_len(number_of_colors)]
  } else {
    # Use rainbow for larger sets
    category_colors <- rainbow(number_of_colors)
  }

  # Name the color vector
  names(category_colors) <- category_values

  # Add color column to metadata
  metadata_df$track_color <- category_colors[metadata_df[[category_to_color_by]]]

  message("Colors generated:")
  for (i in seq_along(category_values)) {
    message("  ", category_values[i], ": ", category_colors[i])
  }

  message("Section 7 complete\n")

} else {
  message("=== SECTION 7: Skipped (COLOR_SCHEMES = FALSE) ===\n")
}
#===============================================================================
# SECTION 8: Create Grouping Columns
#===============================================================================
if (RECOMPUTATION_FLAGS$GROUPING_COLUMNS) {
  message("=== SECTION 8: Create Grouping Columns ===")

  # Validate EXPERIMENTAL_GROUPING_COLUMNS exists
  if (!exists("EXPERIMENTAL_GROUPING_COLUMNS", where = EXPERIMENT_CONFIG)) {
    message("EXPERIMENTAL_GROUPING_COLUMNS not found, skipping grouping columns")
    message("Section 8 complete\n")

  } else {

    grouping_columns_list <- EXPERIMENT_CONFIG$EXPERIMENTAL_GROUPING_COLUMNS

    # Validate it's a named list
    if (!is.list(grouping_columns_list) || is.null(names(grouping_columns_list))) {
      stop("EXPERIMENTAL_GROUPING_COLUMNS must be a named list.\n",
           "Example: list(track_label = c('col1', 'col2'))")
    }

    message("Creating ", length(grouping_columns_list), " grouping columns...")

    # Create each grouping column
    for (group_name in names(grouping_columns_list)) {
      grouping_columns <- grouping_columns_list[[group_name]]

      # Validate grouping columns exist in metadata
      missing_cols <- setdiff(grouping_columns, colnames(metadata_df))
      if (length(missing_cols) > 0) {
        stop("Grouping '", group_name, "' references missing columns:\n",
             "  ", paste(missing_cols, collapse = "\n  "), "\n",
             "Available columns: ", paste(colnames(metadata_df), collapse = ", "))
      }

      # Create grouping column by pasting specified columns
      metadata_df[[group_name]] <- do.call(
        paste,
        c(lapply(metadata_df[, grouping_columns], as.character), sep = "-")
      )

      # Report results
      n_unique <- length(unique(metadata_df[[group_name]]))
      message("  ", group_name, ":")
      message("    Based on: ", paste(grouping_columns, collapse = ", "))
      message("    Unique groups: ", n_unique)
    }

    message("Section 8 complete\n")

  }

} else {
  message("=== SECTION 8: Skipped (GROUPING_COLUMNS = FALSE) ===\n")

}

#===============================================================================
# SECTION 9: Compute Fragment Sizes
#===============================================================================
if (RECOMPUTATION_FLAGS$FRAGMENT_SIZES) {
  message("=== SECTION 9: Compute Fragment Sizes ===")

  # TODO: Implement fragment size calculation
  # This is a placeholder for future implementation

  # Initialize column
  metadata_df$fragment_size <- NA

  message("Fragment size computation not yet implemented")
  message("Section 9 complete\n")

} else {
  message("=== SECTION 9: Skipped (FRAGMENT_SIZES = FALSE) ===\n")
}

#===============================================================================
# SECTION 10: Save Prepared Metadata
#===============================================================================
message("=== SECTION 10: Save Prepared Metadata ===")

# Define output path
output_path <- file.path(
  EXPERIMENT_DIR_PATH,
  "documentation",
  paste0(EXPERIMENT_ID, "_prepared_metadata.csv")
)

# Check if file exists and we're not in dry run
if (file.exists(output_path) && !DRY_RUN) {
  if (interactive()) {
    user_input <- readline(prompt = sprintf(
      "File exists: %s\nOverwrite? (y/n): ",
      basename(output_path)
    ))

    if (tolower(user_input) != "y") {
      message("Save cancelled by user")
      message("Section 10 complete\n")
      stop("Script terminated by user", call. = FALSE)
    }
  } else {
    message("Overwriting existing file: ", basename(output_path))
  }
}

# Write output
if (DRY_RUN) {
  message("[DRY RUN] Would save: ", basename(output_path))
  message("  Rows: ", nrow(metadata_df))
  message("  Columns: ", ncol(metadata_df))
  message("  Column names: ", paste(colnames(metadata_df), collapse = ", "))
} else {
  tryCatch({
    write.csv(metadata_df, output_path, row.names = FALSE)
    message("Saved: ", basename(output_path))
    message("  Rows: ", nrow(metadata_df))
    message("  Columns: ", ncol(metadata_df))

  }, error = function(e) {
    stop("Failed to save prepared metadata: ", e$message)
  })
}

message("Section 10 complete\n")

#===============================================================================
# Metadata Preparation Complete
#===============================================================================
message("=== Metadata Preparation Complete ===")
message("Experiment: ", EXPERIMENT_ID)
message("Output: ", basename(output_path))
message("Ready for downstream analysis and plotting\n")

################################################################################
# End of prepare_experiment_metadata.R
################################################################################
