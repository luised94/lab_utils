###############################################################################
# Load and analyze peak files
################################################################################
# PURPOSE: Load and analyze peak files produced by pipeline test sequences
# Conclusion:
#   = {{fill}}
# USAGE: source("reference_code/pipeline_completion/analyze_peak_files.R")
# DEPENDENCIES: GenomicRanges, rtracklayer
# OUTPUT: {{fill}}
# AUTHOR: LEMR
# DATE: 2025-04-08
# DATE_V1_COMPLETE: {{fill}}
################################################################################
# Bootstrap phase
FUNCTION_FILENAMES <- c("logging", "script_control", "file_operations")
for (function_filename in FUNCTION_FILENAMES) {
    function_filepath <- sprintf("~/lab_utils/core_scripts/functions_for_%s.R", function_filename)
    normalized_path <- normalizePath(function_filepath)
    if (!file.exists(normalized_path)) {
        stop(sprintf("[FATAL] File with functions not found: %s", normalized_path))
    }
    source(normalized_path)
}

# Proceed if packages are installed. Can be disable.
REQUIRED_PACKAGES <- c("rtracklayer", "GenomicRanges", "Gviz")
check_required_packages(REQUIRED_PACKAGES, verbose = TRUE, skip_validation = FALSE)
################################################################################
# Setup genome and feature files
################################################################################
# Ensure supplementary files for plotting genome tracks are present
FILE_GENOME_DIRECTORY <- file.path(Sys.getenv("HOME"), "data", "REFGENS")
FILE_GENOME_PATTERN <- "S288C_refgenome.fna"
FILE_FEATURE_DIRECTORY <- file.path(Sys.getenv("HOME"), "data", "feature_files")
FILE_FEATURE_PATTERN <- "eaton_peaks"
#SGD_FILE_FEATURE_PATTERN <-
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

if (length(REF_GENOME_FILE) == 0) {
    stop(sprintf("No reference genome files found matching pattern '%s' in: %s",
                FILE_GENOME_DIRECTORY,
                FILE_FEATURE_DIRECTORY))
}

if (!file.exists(REF_GENOME_FILE)) {
    stop(sprintf("Reference genome file not accessible: %s", REF_GENOME_FILE[1]))
}

REFERENCE_GENOME_DSS <- Biostrings::readDNAStringSet(REF_GENOME_FILE)

# Create chromosome range
CHROMOSOME_TO_PLOT <- 10
CHROMOSOME_WIDTH <- REFERENCE_GENOME_DSS[CHROMOSOME_TO_PLOT]@ranges@width
CHROMOSOME_ROMAN <- paste0("chr", utils::as.roman(CHROMOSOME_TO_PLOT))

GENOME_RANGE_TO_LOAD <- GenomicRanges::GRanges(
    seqnames = CHROMOSOME_ROMAN,
    ranges = IRanges::IRanges(start = 1, end = CHROMOSOME_WIDTH),
    strand = "*"
)

# Load feature file (annotation)
FEATURE_FILE <- list.files(
    path = FILE_FEATURE_DIRECTORY,
    pattern = FILE_FEATURE_PATTERN,
    full.names = TRUE
)[1]

if (length(FEATURE_FILE) == 0)
{
    warning(sprintf("No feature files found matching pattern '%s' in: %s",
                   FILE_FEATURE_PATTERN,
                   FILE_FEATURE_DIRECTORY))
}

if (!is.null(FEATURE_FILE))
{
    GENOME_FEATURES <- rtracklayer::import(FEATURE_FILE)
    # Convert to chrRoman format
    GenomeInfoDb::seqlevels(GENOME_FEATURES) <- paste0(
        "chr",
        utils::as.roman(gsub("chr", "", GenomeInfoDb::seqlevels(GENOME_FEATURES)))
    )
    # Fitler out the rest of the chromosomes
    # GENOME_FEATURES <- GENOME_FEATURES[seqnames(GENOME_FEATURES) == CHROMOSOME_ROMAN]
    #GENOME_FEATURES <- GenomeInfoDb::keepSeqlevels(GENOME_FEATURES, CHROMOSOME_ROMAN, pruning.mode = "coarse")
}

################################################################################
# Setup directories and file paths
################################################################################

#CONTROL_SAMPLES <- list(
#    "test" = "$HOME/data/241010Bel/alignment/processed_245018_sequence_to_S288C_sorted.bam",
#    "input" = "$HOME/data/241010Bel/alignment/processed_245003_sequence_to_S288C_sorted.bam",
#    "reference" = "$HOME/data/100303Bel/alignment/processed_034475_sequence_to_S288C_sorted.bam"
#)

DATA_DIRECTORY <- file.path(Sys.getenv("HOME"), "data", "preprocessing_test")
SUPPORTED_FILE_TYPES <- c(
  FASTQ = "fastq",
  BAM = "bam",
  BW = "bw",
  NARROW_PEAK = "narrowPeak",
  BED = "bed",
  XLS = "xls",
  GAPPED_PEAK = "gappedPeak",
  BROAD_PEAK = "broadPeak"
)

# Preallocate list with file paths
file_paths_by_type <- setNames(vector("list", length(SUPPORTED_FILE_TYPES)), names(SUPPORTED_FILE_TYPES))

for (file_type in names(SUPPORTED_FILE_TYPES)) {
    extension <- SUPPORTED_FILE_TYPES[[file_type]]
    file_pattern <- paste0(".", extension, "$")
    message(sprintf("Searching for '%s' files (pattern: %s)...", file_type, file_pattern))

    found_files <- list.files(
        path = DATA_DIRECTORY,
        pattern = file_pattern,
        recursive = TRUE,
        full.names = TRUE,
        include.dirs = FALSE
    )

    if (length(found_files) == 0) {
        warning(sprintf("No files found for type: %s", file_type))
        next
    }

    message(sprintf("Found %s files for '%s' files ...", length(found_files), file_type))
    file_paths_by_type[[file_type]] <- found_files
}

# Remove empty entries if warnings are acceptable
# If no files are found, next assertion will stop
file_paths_by_type <- Filter(length, file_paths_by_type)
if (length(file_paths_by_type) == 0) {
    stop("No valid files found in directory: ", DATA_DIRECTORY)
}

################################################################################
# Load and process peak files
################################################################################
# Load PEAK_FILE_COLUMNS variable
source(normalizePath("~/lab_utils/core_scripts/peak_file_columns_by_type.R"))
stopifnot("Problem loading PEAK_FILE_COLUMNS variable from R file." = exists("PEAK_FILE_COLUMNS"))

# Setup
#PEAK_EXTENSIONS <- c("NARROW_PEAK", "BED", "XLS", "GAPPED_PEAK", "BROAD_PEAK")
#PEAK_FILES <- file_paths_by_type[PEAK_EXTENSIONS]
EXTENSION_TO_LOAD <- "XLS"
EXTENSIONS_TO_REMOVE <- paste0(".", tolower(EXTENSION_TO_LOAD))
ESCAPED_EXTENSIONS <- gsub(".", "\\.", EXTENSIONS_TO_REMOVE, fixed = TRUE)
XLS_FILES <- file_paths_by_type[[EXTENSION_TO_LOAD]]

if (length(XLS_FILES) == 0) {
    stop("No xls peak files found.")
}

#-------------------------------------------------------------------------------
# EXTRACT METADATA FROM FILENAMES
#-------------------------------------------------------------------------------
# Extract base filenames without path and extensions
file_basenames <- basename(XLS_FILES)
filenames_without_extension <- gsub(ESCAPED_EXTENSIONS, "", file_basenames)
clean_filenames <- gsub("_peaks$", "", filenames_without_extension) # Remove _peaks suffix added by MACS2

# Split filenames by underscores to extract metadata components
filename_parts <- strsplit(clean_filenames, split = "_")
max_parts <- max(sapply(filename_parts, length))

# Pad each split filename to ensure uniform length for data frame creation
padded_filename_parts <- lapply(filename_parts, function(parts) {
  padding_needed <- max_parts - length(parts)
  if (padding_needed > 0) {
    return(c(parts, rep("none", padding_needed)))
  } else {
    return(parts)
  }
})

# Convert to data frame
metadata_df <- as.data.frame(do.call(rbind, padded_filename_parts))
ORIGINAL_COLUMN_COUNT <- ncol(metadata_df)

# Count "none" values in each row (vectorized approach)
metadata_df$none_count <- apply(metadata_df, 1, function(row) {
  sum(grepl("none", row))
})

# Determine BAM type based on position in filename
ADDITIONAL_NUMBER_OF_METADATA_COLUMNS <- 2
metadata_df$bam_type <- sapply(1:nrow(metadata_df), function(row_idx) {
  # Calculate position of BAM type in the row
  bam_type_col_idx <- ORIGINAL_COLUMN_COUNT - ADDITIONAL_NUMBER_OF_METADATA_COLUMNS - metadata_df$none_count[row_idx]
  if (bam_type_col_idx > 0 && bam_type_col_idx <= ORIGINAL_COLUMN_COUNT) {
    return(as.character(metadata_df[row_idx, bam_type_col_idx]))
  } else {
    return("unknown")  # Fallback if position is invalid
  }
})

# Extract the last two non-"none" columns as input_type and peak_type
input_and_peak_cols <- data.frame(
  input_type = character(nrow(metadata_df)),
  peak_type = character(nrow(metadata_df))
)

for (row_idx in 1:nrow(metadata_df)) {
  end_idx <- ORIGINAL_COLUMN_COUNT - metadata_df$none_count[row_idx]
  if (end_idx >= 2) {
    start_idx <- end_idx - 1
    input_and_peak_cols[row_idx, ] <- metadata_df[row_idx, start_idx:end_idx]
  }
}

# COMBINE ALL METADATA INTO FINAL DATA FRAME
# Assuming first two columns are sample ID and condition
final_metadata_df <- data.frame(
  sample_type = metadata_df[, 1],
  condition_idx = metadata_df[, 2],
  bam_type = metadata_df$bam_type,
  input_type = input_and_peak_cols$input_type,
  peak_type = input_and_peak_cols$peak_type,
  file_path = XLS_FILES,
  stringsAsFactors = FALSE
)

# Assertions
# TODO: Need to see if I should extract the number into variables to improve readability. //
stopifnot(
    "No NAs found in the final_metadata_df." = !any(is.na(final_metadata_df)),
    "Metadata has expected dimensions." =
        nrow(final_metadata_df) == length(XLS_FILES) &&
        ncol(final_metadata_df) == ADDITIONAL_NUMBER_OF_METADATA_COLUMNS + 2 + 1 + 1,
    "No duplicate samples." = !any(duplicated(final_metadata_df))
)

# Determine unique values for each categories
message("Determining unique categories...")
metadata_columns_vec <- setdiff(names(final_metadata_df), "file_path")
unique_categories_lst <- lapply(metadata_columns_vec, function(col_name) {
  unique(final_metadata_df[[col_name]])
})
names(unique_categories_lst) <- metadata_columns_vec
lapply(metadata_columns_vec, function(col_name) {
  message(sprintf("--- Column %s ---\n  %s",
            col_name,
            paste(unique_categories_lst[[col_name]], collapse = ",")))
})

# Define the levels of the columns that should be turned into factors
# See the output of the previous lapply statement to manually define
# the order for the factors. This will control the plotting order.
factor_levels_list <- list(
  bam_type = c("raw", "deduped", "shifted", "blFiltered"),
  peak_type = c("narrow", "broad")
)
stopifnot(all(names(factor_levels_list) %in% names(final_metadata_df)))
for (column_name in names(factor_levels_list)){
  expected_levels <- factor_levels_list[[column_name]]
  actual_levels <- unique(final_metadata_df[[column_name]])
  missing_levels <- setdiff(expected_levels, actual_levels)
  if (length(missing_levels) > 0) {
    stop(sprintf(
        "Validation failed for column '%s':\n  Levels defined but missing: %s\n  Available levels: %s",
        column_name,
        paste(missing_levels, collapse = ", "),
        paste(actual_levels, collapse = ", ")
        ),
        call. = FALSE
    )
  }
} # end validation for loop
for (col in names(factor_levels_list)) {
  final_metadata_df[[col]] <- factor(
    x = final_metadata_df[[col]],
    levels = factor_levels_list[[col]],
    ordered = TRUE
  )
}
message("User specified columns converted to factors...")

##################################################################################
# MAIN
##################################################################################
# Create list of unique categories automatically
#UNIQUE_METADATA_CATEGORIES <- lapply(colnames(final_metadata), function(col_name) {
#  unique(final_metadata[[col_name]])
#})
#
## Name the list elements using the column names
#names(UNIQUE_METADATA_CATEGORIES) <- colnames(final_metadata)

# Create output directory if it doesn't exist
PLOT_OUTPUT_DIR <- file.path(Sys.getenv("HOME"), "data", "preprocessing_test", "plots")
dir.create(PLOT_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Configuration ---
MIN_ROW_COUNT <- 1
MAX_ROW_COUNT <- nrow(final_metadata_df)
ANALYSIS_COLUMNS <- c("number_of_peaks", "peak_width_distribution", "overlap_with_eaton")

DEBUG_MODE <- TRUE
# Prevent time-consuming code from running unless it has to
if (interactive() && !exists("summary_statistics_df")){
cat("Running the summary statistics computation...\n")
# Code block that is "cached" --------
# 1. Define types and order
SUMMARY_COLS <- c("sample_id", "file_path", "num_peaks",
                 "percent_recovered", "percent_enriched",
                 "width_mean", "width_median")

# 2. Initialize with correct types
summary_statistics_df <- data.frame(
  sample_id = character(MAX_ROW_COUNT),
  file_path = character(MAX_ROW_COUNT),
  num_peaks = integer(MAX_ROW_COUNT),
  percent_recovered = numeric(MAX_ROW_COUNT),
  percent_enriched = numeric(MAX_ROW_COUNT),
  width_mean = numeric(MAX_ROW_COUNT),
  width_median = numeric(MAX_ROW_COUNT),
  stringsAsFactors = FALSE
)[, SUMMARY_COLS]  # Enforce order

ESTIMATED_CHROMOSOMES_PER_SAMPLE <- length(names(REFERENCE_GENOME_DSS))
chromosome_distribution_list <- vector("list", MAX_ROW_COUNT)
peak_statistics_list <- vector("list", MAX_ROW_COUNT)

message(sprintf("Starting processing for %d records...", MAX_ROW_COUNT))
for (row_index in seq_len(MAX_ROW_COUNT)) {
  # Pre-loop validation
  stopifnot(
    row_index >= MIN_ROW_COUNT,
    row_index <= MAX_ROW_COUNT
  )
  # Add a small separator for clarity in logs
  message("====================")
  message(sprintf("Processing row %d/%d", row_index, MAX_ROW_COUNT))

  # Extract the current row's data
  current_metadata <- final_metadata_df[row_index, ]
  if (DEBUG_MODE) {
    message("Row Details:")
    invisible(lapply(names(current_metadata), function(column_name) {
        column_data <- current_metadata[[column_name]]
       message(sprintf("| %s: %s - %s", column_name, column_data, typeof(column_data)))
    }))
  }

  xls_file_path <- as.character(current_metadata$file_path)
  current_sample_id <- paste0(current_metadata$sample_type, as.character(current_metadata$condition_idx))

  # --- Load xls file data ---
  message("--- Loading xls data ---")
  stopifnot(
    "Expect only one file path in xls_file_path every iteration." = length(xls_file_path) == 1,
    "xls_file_path variable is empty." = nzchar(xls_file_path),
    "xls_file_path variable for current metadata does not exist." = file.exists(xls_file_path)
  )

  file_lines <- readLines(xls_file_path)
  # Note: There is one line in the header that is empty. //
  # This should account for it but there may be issues.  //
  comment_lines <- grep("^#", file_lines)
  # Determine the header line (first non-comment line)
  if (length(comment_lines) > 0) {
    header_line <- max(comment_lines) + 1
  } else {
    header_line <- 1  # No comments, header is first line
  }
  message(sprintf("  Current sample id: %s", current_sample_id))
  message(sprintf("  Number of comment lines: %s", length(comment_lines)))
  message(sprintf("  Header line identified: %s", header_line))

  # Check if file has content after comments
  if (header_line >= length(file_lines)) {
    message(paste("  File has no data after comments:", xls_file_path))
    warning(paste("File has no data after comments:", xls_file_path))
    message("  Creating placeholder data...")
    new_row <- data.frame(
      sample_id = current_sample_id,
      file_path = xls_file_path,
      num_peaks = 0L,
      percent_recovered = NA_real_,
      percent_enriched = NA_real_,
      width_mean = NA_real_,
      width_median = NA_real_,
      stringsAsFactors = FALSE
    )
    stopifnot(all(names(summary_statistics_df) %in% names(new_row)))
    new_row <- new_row[, names(summary_statistics_df)]
    stopifnot(identical(names(new_row), names(summary_statistics_df)))
    summary_statistics_df[row_index, names(new_row)] <- new_row

    # Create minimal placeholder data for distributions
    chromosome_distribution_list[[row_index]] <- data.frame(
      sample_id = current_sample_id,
      file_path = xls_file_path,
      chromosome = NA_character_,
      count = NA_integer_,
      stringsAsFactors = FALSE
    )
    peak_statistics_list[[row_index]] <- data.frame(
      sample_id = current_sample_id,
      file_path = xls_file_path,
      peak_id = NA_character_,
      fold_enrichment = NA_real_,
      qvalue = NA_real_,
      width = NA_real_,
      stringsAsFactors = FALSE
    )
    message("  Skipping to next iteration...")
    next  # Skip to next iteration
  }

  # Get the header and data lines
  header <- file_lines[header_line]
  data_lines <- file_lines[(header_line+1):length(file_lines)]

  # Check if there's any data
  if (length(data_lines) == 0) {
    warning(paste("File has header but no data:", file_path))
    next
  }

  # --- Extract statistics from xls file ---
  message("--- Extracting statistics from xls file ---")
  xls_peak_df <- read.delim(text = data_lines, header = FALSE, sep = "\t")
  message(sprintf(
    "  Loaded %s peaks (cols: %s) from %s",
    nrow(xls_peak_df),
    ncol(xls_peak_df),
    basename(xls_file_path)
  ))

  file_ext <- tools::file_ext(xls_file_path)
  format_key <- toupper(file_ext)
  # --- Handle --call-summits command: Call summit adds one more column to the output ---
  # Solution: Use a different set of columns names
  if (any(grepl("--call-summits",
                file_lines[1:max(comment_lines)],
                ignore.case = TRUE))) {

    format_key <- paste0(format_key, "_SUMMITS")
  }
  SUPPORTED_FORMATS <- paste(tolower(names(PEAK_FILE_COLUMNS)), collapse=", ")
  # Check if we have a definition for this file type
  if (!format_key %in% names(PEAK_FILE_COLUMNS)) {
    stop(paste0("Unsupported file format: ", file_ext,
                ". Supported formats are: ",
                SUPPORTED_FORMATS))
  }
  col_names <- PEAK_FILE_COLUMNS[[format_key]]
  names(xls_peak_df) <- col_names
  # --- Error handling chunk: Ensure loaded dataframe has ---
  # --- expected properties for genomic range data. ----
  has_chr_prefix <- grepl("^chr", xls_peak_df$chromosome)
  stopifnot(
    "xls_peak_df is empty." = nrow(xls_peak_df) > 0,
    "Missing chromosome column in xls_peak_df." = "chromosome" %in% names(xls_peak_df),
    "Missing start position column in xls_peak_df." = "start" %in% names(xls_peak_df),
    "Missing end position column in xls_peak_df." = "end" %in% names(xls_peak_df),
    "Some chromosome values in xls_peak_df do not follow chr<roman_num> format." = all(has_chr_prefix),
    "xls_peak_df has negative start positions" = all(xls_peak_df$start >= 0),
    "xls_peak_df has rows where end < start" = all(xls_peak_df$end > xls_peak_df$start)
  )

  lapply(names(xls_peak_df), function(column_name){
    if (any(is.na(xls_peak_df[[column_name]]))){
      warning(sprintf("Some missing values in %s column.", column_name))
    }
  })
  # -----------

  message(sprintf(" Error handling chunk passed for: %s", xls_file_path))
  message(sprintf("  Peak dataframe details from file path: %s", xls_file_path))
  message("~~~~~~~~~~~~~~~~~~~~")
  print(head(xls_peak_df))
  message("~~~~~~~~~~~~~~~~~~~~")

  # --- Create GenomicRanges ---
  # Create GRanges object with core coordinates
  gr <- GenomicRanges::makeGRangesFromDataFrame(xls_peak_df, keep.extra.columns = TRUE)
  message("  Assigned Genomic Ranges...")

  # --- Bed File Analysis Basic Statistics ---
  # Validate inputs
  if (length(gr) == 0) {
    warning("No peaks found in the current dataset")
    next
  }

  if (length(GENOME_FEATURES) == 0) {
    message("[ERROR] Reference genome features object is empty.")
    stop("Reference genome is empty and required for downstream analysis.")
  }

  message("  Calculating overlaps...")
  message("  Number of called peaks: ", length(gr))
  message("  Number of reference features: ", length(GENOME_FEATURES))

  # --- Peak widths: What is the width of the called peaks? ---
  # Is it consistent with previous experiments and DNA binding behavior?
  peak_widths <- xls_peak_df$end - xls_peak_df$start
  # --- Chromosome distribution: How many peaks in each chromosome? ---
  # chrom_counts <- as.data.frame(table(factor(GenomicRanges::seqnames(gr), levels = GenomicInfoDb::seqlevels(gr)))
  chrom_counts <- as.data.frame(table(
    # Requires chromosomes of the reference genome to be in order
    factor(xls_peak_df$chromosome, levels=names(REFERENCE_GENOME_DSS))
  ))

  # --- Enrichment: What percent of called peaks overlap with reference features ---
  # Is it consistent with previous experiments of ORC and MCM?
  # Calculate overlaps (logical vector: TRUE = overlap exists)
  # Use overlapsAny to avoid double counting.
  current_enrichment_overlaps <- IRanges::overlapsAny(gr, GENOME_FEATURES)
  # Count overlapping peaks
  current_overlapping_peaks <- sum(current_enrichment_overlaps)
  # Percent of called peaks overlapping reference features
  current_percent_enriched <- (current_overlapping_peaks / length(gr)) * 100

  # --- Recovery: What percent of reference features are recovered by called peaks ---
  # Calculate overlaps (logical vector: TRUE = overlap exists)
  current_recovery_overlaps <- IRanges::overlapsAny(GENOME_FEATURES, gr)
  # Count recovered reference features
  current_recovered_features <- sum(current_recovery_overlaps)
  # Percent of reference features overlapped by called peaks
  current_percent_recovered <- (current_recovered_features / length(GENOME_FEATURES)) * 100

  # --- Validation check: Percent should not exceed 100 ---
  if (current_percent_enriched > 100) {
    warning_message <- sprintf("Enrichment percentage exceeds 100% for:\n  [ %s ]",
                               current_sample_id)
    warning(warning_message)
  }
  if (current_percent_recovered > 100) {
    warning_message <- sprintf("Recovery percentage exceeds 100% for:\n  [ %s ]",
                               current_sample_id)
    warning(warning_message)
  }

  stopifnot(
    "Percent recovered exceeds 100." = current_percent_recovered <= 100,
    "Percent enriched exceeds 100." = current_percent_enriched <= 100
  )

  message("  Calculated statistics based on dataframe...")
  # --- Detailed Diagnostics ---
  message("  Enrichment: ", current_overlapping_peaks, "/", length(gr),
          " called peaks overlap with reference features (",
          formatC(current_percent_enriched, digits=2, format="f"), "%)")
  message("  Recovery: ", current_recovered_features, "/", length(GENOME_FEATURES),
          " reference features recovered by called peaks (",
          formatC(current_percent_recovered, digits=2, format="f"), "%)")

  stopifnot(
    "Percent recovered exceeds 100." = current_percent_recovered <= 100,
    "Percent enriched exceeds 100." = current_percent_enriched <= 100
  )
  # --- Store Results ---
  number_of_peaks <- nrow(xls_peak_df)
  # Critical: Order of assignment in the data.frame call has to be the same!!!
  new_row <- data.frame(
    sample_id = current_sample_id,
    file_path = xls_file_path,
    num_peaks = number_of_peaks,
    percent_recovered = current_percent_recovered,
    percent_enriched = current_percent_enriched,
    width_mean = mean(peak_widths),
    width_median = median(peak_widths),
    stringsAsFactors = FALSE
  )
  stopifnot(all(names(summary_statistics_df) %in% names(new_row)))
  new_row <- new_row[, names(summary_statistics_df)]
  stopifnot(identical(names(new_row), names(summary_statistics_df)))
  summary_statistics_df[row_index, names(new_row)] <- new_row

  message("  Added statistics to summary df...")

  # Create minimal placeholder data for distributions
  chromosome_distribution_list[[row_index]] <- data.frame(
    sample_id = current_sample_id,
    file_path = xls_file_path,
    chromosome = chrom_counts$Var1,
    count = chrom_counts$Freq,
    stringsAsFactors = FALSE
  )
  message("  Added statistics to chromosome distribution list...")
  peak_statistics_list[[row_index]] <- data.frame(
    sample_id = rep(current_sample_id, number_of_peaks),
    file_path = rep(xls_file_path, number_of_peaks),
    peak_id = seq_len(number_of_peaks),
    xls_peak_df[c("fold_enrichment", "qvalue")],
    width = peak_widths,
    stringsAsFactors = FALSE
  )
  message("  Added statistics to peak statistics list...")

  message("  Finishing loop iteration...")
  message("====================")
} # End of for loop
message("Processing finished.")
# Convert lists to DF
peak_statistics_df <- do.call(rbind, peak_statistics_list)
chromosome_distribution_df <- do.call(rbind, chromosome_distribution_list)

# Remove any rows where all values are NA.
# Peaks will be zero in some cases which shouldnt cause drop.
chromosome_distribution_df <- chromosome_distribution_df[!apply(is.na(chromosome_distribution_df), 1, all), ]
peak_statistics_df <- peak_statistics_df[!apply(is.na(peak_statistics_df), 1, all), ]
summary_statistics_df <- summary_statistics_df[!apply(is.na(summary_statistics_df), 1, all), ]

# Add metadata to distribution DFs for easier plotting
#METADATA_COLS_TO_KEEP <- setdiff(names(final_metadata_df), "file_paths")

# Explicit merges
chromosome_distribution_df <- merge(
  chromosome_distribution_df,
  final_metadata_df,
  by = "file_path"
)
peak_statistics_df <- merge(
  peak_statistics_df,
  final_metadata_df,
  by = "file_path"
)
summary_statistics_df <- merge(
  summary_statistics_df,
  final_metadata_df,
  by = "file_path"
)
message("Merge complete...")
#stopifnot(
  #"Percent recovered exceeds 100." = all(summary_statistics_df$percent_recovered <= 100),
  #"Percent enriched exceeds 100." = all(summary_statistics_df$percent_enriched <= 100)
#)

# End of "cached" code ----------
cat("Computation complete! -- summary_statistics_df is now available.\n")

# Cache the for loop ###########
} else if (interactive()) {
  cat("Skipping summary statistics computation: 'summary_statistics_df' already exists.\n",
      "To force a rerun, remove the variable with: rm(summary_statistics_df)\n")
}

# --- Plot results ---
# You can modify and copy to the interactive session in the cluster
# to avoid triggering git. Otherwise I would have to create a config file.
# Similar to how I control the r scripts in the core_scripts dir.
# TODO: May create a config in the future but this is just a test set of scripts
library(tidyverse)
library(ggplot2)

# Configure sample selection:
all_sample_ids_chr <- unique(summary_statistics_df$sample_id)
# Subset, grep or all as examples
sample_ids_to_plot_chr <- all_sample_ids_chr
# sample_ids_to_plot <- grep("input", all_sample_ids, value = TRUE)
sample_ids_to_plot_chr <- all_sample_ids_chr[2]

# Reference values for the plots
reference_peak_count_int <- length(GENOME_FEATURES)
recovery_reference_percent <- 100

message("====================")
for (current_sample_id in sample_ids_to_plot_chr) {
  message("--- Plot sample id for loop ---")
  is_row_with_sample_id_bool <- summary_statistics_df$sample_id == current_sample_id
  current_sample_subset_df <- summary_statistics_df[is_row_with_sample_id_bool, ]

  # Skip if no rows in the subset dataframe
  if (nrow(current_sample_subset_df) == 0) {
    warning(sprintf("Skipping empty sample: %s", current_sample_id))
    next

  }

  message(sprintf(
    paste0("  Subset dataframe with sample_id: %s\n",
           "  Number of rows in subset df: %s"),
    current_sample_id,
    nrow(current_sample_subset_df)
  ))

  # Create processing group factor
  current_sample_subset_df$processing_group <- factor(
    paste(current_sample_subset_df$bam_type, current_sample_subset_df$peak_type, sep = " + ")
  )

  # Sort processing_group by median recovery % (descending)
  # Create sorted version of processing_group for recovery
  current_sample_subset_df <- current_sample_subset_df %>%
    mutate(
      processing_group_srt_rec = fct_reorder(
        processing_group,
        percent_recovered,
        .fun = median,
        .desc = TRUE
      ),
      processing_group_srt_enr = fct_reorder(
        processing_group,
        percent_enriched,
        .fun = median,
        .desc = TRUE
      )
    )
  # --- How many peaks were recovered ---
  peak_recovery_plot <- ggplot(
    current_sample_subset_df,
    aes(x = processing_group_srt_rec, y = percent_recovered,
        color = input_type, group = input_type)
    ) +
    geom_hline(yintercept = recovery_reference_percent,
               linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_text(
      aes(x = Inf, y = recovery_reference_percent,
          label = paste("Reference =", recovery_reference_percent, "%")),
      vjust = -0.5, hjust = 1.1,
      color = "black", size = 3.5) +
    geom_line(linewidth = 0.7, alpha = 0.5) +
    geom_point(size = 3) +
    geom_text(
      aes(label = round(percent_recovered, 1)),
          vjust = -1, size = 3, show.legend = FALSE) +
      labs(title = "Peak Recovery by Processing Method",
           subtitle = paste("Sample:", current_sample_id),
           x = "Processing Method (sorted by recovery %)",
           y = "Peaks Recovered (%)",
           color = "Input Control") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # --- How many peaks are enriched  ---
  peak_enrichment_plot <- ggplot(
    current_sample_subset_df,
    aes(x = processing_group_srt_enr, y = percent_enriched,
        color = input_type, group = input_type)
    ) +
    geom_hline(yintercept = recovery_reference_percent,
               linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_text(
      aes(x = Inf, y = recovery_reference_percent,
          label = paste("Reference =", recovery_reference_percent, "%")),
      vjust = -0.5, hjust = 1.1,
      color = "black", size = 3.5) +
    geom_line(linewidth = 0.7, alpha = 0.5) +
    geom_point(size = 3) +
    geom_text(
      aes(label = round(percent_enriched, 1)),
          vjust = -1, size = 3, show.legend = FALSE) +
      labs(title = "Peak Enrichment by Processing Method",
           subtitle = paste("Sample:", current_sample_id),
           x = "Processing Method (sorted by recovery %)",
           y = "Peaks Enriched (%)",
           color = "Input Control") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # --- How many peaks were called ---
  peak_count_plot <- ggplot(current_sample_subset_df,
    aes(x = processing_group, y = num_peaks, fill = input_type)) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_hline(yintercept = reference_peak_count_int,
               linetype = "dashed",
               linewidth = 0.8,
               color = "black") +
    geom_text(
      aes(x = Inf, y = reference_peak_count_int,
          label = paste("Eaton Peaks =", reference_peak_count_int)),
      vjust = -0.5, hjust = 1.1,
      color = "black", size = 3.7, fontface = "bold") +
    labs(title = "Number of Peaks Called",
         subtitle = paste("Sample:", current_sample_id),
         x = "Processing Method", y = "Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  #peak_width_stats_plot <- ggplot(current_sample_subset_df) +
  #  geom_point(aes(x = processing_group, y = width_mean,
  #                 fill = input_type, color = "Mean"), size = 3) +
  #  geom_point(aes(x = processing_group, y = width_median,
  #                 fill = input_type, color = "Median"), size = 3) +
  #  labs(
  #    title = "Peak Width Statistics\n",
  #    subtitle = paste("Sample:", current_sample_id),
  #    x = "Processing Method + Peak Type",
  #    y = "Width (bp)",
  #    color = "Statistic") +
  #  theme_minimal() +
  #  theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Get all ggplot objects from environment
  plot_object_names <- ls(pattern = "_plot$", envir = .GlobalEnv)
  for (current_plot_name in plot_object_names) {
    message("  --- Plotting variables ---")
    message(sprintf("  Plotting: %s", current_plot_name))
    current_plot_object <- get(current_plot_name, envir = .GlobalEnv)
    if (!inherits(current_plot_object, "ggplot")) {
      warning("Skipping ", current_plot_name, " - not a ggplot object")
      next
    }

    plot_output_path <- file.path(
      PLOT_OUTPUT_DIR,
      paste0(current_plot_name, "_", current_sample_id, ".svg")
    )
    message(sprintf("  Saving to: %s", plot_output_path))
    #svglite::svglite(
    #    filename = plot_output_path,
    #    width = 10,
    #    height = 8,
    #    bg = "white"
    #)
    print(current_plot_object)
    #dev.off()
    readline(prompt = "Press [enter] to continue plot")
  }

  # Optional: Add interactive pause between samples
  # invisible(readline(prompt = "Press [enter] to continue to other sample"))
}
message("====================")
