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
final_metadata <- data.frame(
  sample_type = metadata_df[, 1],
  condition_idx = metadata_df[, 2],
  bam_type = metadata_df$bam_type,
  input_type = input_and_peak_cols$input_type,
  peak_type = input_and_peak_cols$peak_type,
  file_path = XLS_FILES,
  stringsAsFactors = FALSE
)

# Assertions
#stopifnot(
#    "No samples are NA." = !any(is.na(final_metadata$sample)),
#    "Metadata has expected dimensions." =
#        nrow(final_metadata) == NUMBER_OF_FILES &&
#        ncol(final_metadata) == NUMBER_OF_COLUMNS,
#    "No duplicate samples." = !any(duplicated(final_metadata))
#)

# Add factor levels to metadata frame and define unique categories
# Convert all columns to factors
#final_metadata[] <- lapply(final_metadata, as.factor)

# Define factor levels for each column (customize as needed)
#COLUMN_LEVELS <- list(
#  sample_group = c("input", "reference", "test"),
#  alignment_processing = c("raw", "deduped", "shifted"),
#  peak_calling_mode = c("narrow", "broad", "auto"),
#  significance_metric = c("p", "q"),
#  input_control_type = c("noInput", "withInput"),
#  output_category = "peaks",  # Single level since constant
#  peak_detail_level = c("peaks", "summits")
#)

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
MAX_ROW_COUNT <- nrow(final_metadata)
ANALYSIS_COLUMNS <- c("number_of_peaks", "peak_width_distribution", "overlap_with_eaton")

DEBUG_MODE <- TRUE
summary_statistics_df <- data.frame(
  sample_id = character(MAX_ROW_COUNT),
  file_path = character(MAX_ROW_COUNT),
  num_peaks = integer(MAX_ROW_COUNT),
  percent_recovered = numeric(MAX_ROW_COUNT),
  percent_enriched = numeric(MAX_ROW_COUNT),
  width_mean = numeric(MAX_ROW_COUNT),
  width_median = numeric(MAX_ROW_COUNT),
  stringsAsFactors = FALSE
)

ESTIMATED_CHROMOSOMES_PER_SAMPLE <- length(names(REFERENCE_GENOME_DSS))
chromosome_distribution_list <- vector("list", MAX_ROW_COUNT)
peak_width_list <- vector("list", MAX_ROW_COUNT)

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
  current_metadata <- final_metadata[row_index, ]
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
  if (header_line > length(file_lines)) {
    message(paste("  File has no data after comments:", xls_file_path))
    warning(paste("File has no data after comments:", xls_file_path))
    message("  Creating placeholder data...")
    summary_statistics_df[row_index, ] <- data.frame(
      sample_id = current_sample_id,
      file_path = xls_file_path,
      num_peaks = 0L,
      width_mean = NA_real_,
      width_median = NA_real_,
      percent_recovered = NA_real_,
      percent_enriched = NA_real_,
      stringsAsFactors = FALSE
    )
    # Create minimal placeholder data for distributions
    chromosome_distribution_list[[row_index]] <- data.frame(
      sample_id = current_sample_id,
      file_path = xls_file_path,
      chromosome = NA_character_,
      count = NA_integer_,
      stringsAsFactors = FALSE
    )
    peak_width_list[[row_index]] <- data.frame(
      sample_id = current_sample_id,
      file_path = xls_file_path,
      peak_id = NA_character_,
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
  # Handle --call-summits command
  # Finds the specific instance in the command line of the header
  # \\b matches whole word.
  if (any(grepl("^# Command:.*--call-summits\\b",
                file_lines[1:header_start],
                ignore.case = TRUE))) {

    format_key <- paste0(format_key, "_SUMMITS")
  }
  SUPPORTED_FORMATS <- paste(tolower(names(PEAK_FILE_COLUMNS)), collapse=", ")
  # Check if we have a definition for this file type
  f (!format_key %in% names(PEAK_FILE_COLUMNS)) {
    stop(paste0("Unsupported file format: ", file_ext,
                ". Supported formats are: ",
                SUPPORTED_FORMATS))
  }
  col_names <- PEAK_FILE_COLUMNS[[format_key]]
  names(xls_peak_df) <- col_names
  # Error handling chunk ----
  has_chr_prefix <- grepl("^chr", xls_peak_df$chromosome)
  stopifnot(
    "xls_peak_df is empty." = nrow(xls_peak_df) > 0,
    "Missing chromosome column in xls_peak_df." = "chromosome" %in% names(xls_peak_df),
    "Missing start position column in xls_peak_df." = "start" %in% names(xls_peak_df),
    "Missing end position column in xls_peak_df." = "end" %in% names(xls_peak_df),
    "Some chromosome values in xls_peak_df do not follow chr<roman_num> format." = any(!has_chr_prefix),
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

  ## --- Create GenomicRanges ---
  ## Create GRanges object with core coordinates
  #gr <- GenomicRanges::GRanges(
  #  seqnames = xls_peak_df$chromosome,
  #  ranges = IRanges::IRanges(
  #    # Convert from 0-based to 1-based
  #    start = xls_peak_df$start + 1L,
  #    end = xls_peak_df$end
  #  ),
  #  strand = "*"
  #  #strand = if("strand" %in% names(xls_peak_df)) xls_peak_df$strand else "*"
  #)

  ## Add metadata columns (everything except chr, start, end)
  #meta_cols <- setdiff(names(xls_peak_df), c("chromosome", "start", "end", "strand"))
  #if (length(meta_cols) > 0) {
  #  #message("Metadata columns:\n", paste(meta_cols, collapse=", "))
  #  GenomicRanges::mcols(gr) <- xls_peak_df[, meta_cols, drop = FALSE]
  #}

  ## --- Bed File Analysis Basic Statistics ---
  #peak_widths <- xls_peak_df$end - xls_peak_df$start
  #chrom_counts <- as.data.frame(table(
  #  # Requires chromosomes of the reference genome to be in order
  #  factor(xls_peak_df$chromosome, levels=names(REFERENCE_GENOME_DSS))
  #))
  ## Calculate overlaps (logical vector: TRUE = overlap exists)
  #overlaps_logical <- IRanges::overlapsAny(gr, GENOME_FEATURES)
  ## Percent of gr1 ranges overlapping gr2
  #percent_enriched <- ( sum(overlaps_logical) / length(gr) ) * 100

  ## Calculate overlaps (logical vector: TRUE = overlap exists)
  #overlaps_logical <- IRanges::overlapsAny(GENOME_FEATURES, gr)
  #percent_recovered <- ( sum(overlaps_logical) / length(GENOME_FEATURES) ) * 100

  ## --- Store Results ---
  #summary_statistics_df[row_index, ] <- data.frame(
  #  sample_id = current_sample_id,
  #   file_path = xls_file_path,
  #  num_peaks = nrow(xls_peak_df),
  #  width_mean = mean(peak_widths),
  #  width_median = median(peak_widths),
  #  percent_recovered = percent_recovered,
  #  percent_enriched = percent_enriched,
  #  stringsAsFactors = FALSE
  #)
  ## Create minimal placeholder data for distributions
  #chromosome_distribution_list[[row_index]] <- data.frame(
  #  sample_id = current_sample_id,
  #    file_path = xls_file_path,
  #  chromosome = chrom_counts$Var1,
  #  count = chrom_counts$Freq,
  #  stringsAsFactors = FALSE
  #)
  #peak_width_list[[row_index]] <- data.frame(
  #  sample_id = current_sample_id,
  #    file_path = xls_file_path,
  #  peak_id = seq_len(length(peak_widths)),
  #  width = peak_widths,
  #  stringsAsFactors = FALSE
  #)

  message("====================")
} # End of for loop
message("Processing finished.")
stop("Breakpoint...")
# Convert lists to DF
peak_width_distribution_df <- do.call(rbind, peak_width_list)
chromosome_distribution_df <- do.call(rbind, chromosome_distribution_list)

# Optional: Trim unused preallocated rows
chromosome_distribution_df <- chromosome_distribution_df[!is.na(chromosome_distribution_df$sample_id), ]
peak_width_distribution_df <- peak_width_distribution_df[!is.na(peak_width_distribution_df$sample_id), ]
summary_statistics_df <- summary_statistics_df[!is.na(summary_statistics_df$sample_id), ]

# Add metadata to distribution DFs for easier plotting
#METADATA_COLS_TO_KEEP <- setdiff(names(final_metadata), "file_paths")

# Explicit merges
chromosome_distribution_df <- merge(
  chromosome_distribution_df,
  final_metadata,
  by = "file_path"
)
peak_width_distribution_df <- merge(
  peak_width_distribution_df,
  final_metadata,
  by = "file_path"
)
summary_statistics_df <- merge(
  summary_statistics_df,
  final_metadata,
  by = "file_path"
)

# --- Plot results ---
