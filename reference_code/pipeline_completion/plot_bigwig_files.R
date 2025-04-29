###############################################################################
# Plot bigwig files
################################################################################
# PURPOSE: Plot bigwig files generated during the pipeline completion tests
# Conclusion: 
#   = The plots look the same for the counts and the post processing. Unless I made a mistake making the bigwig files.
# USAGE: source("reference_code/pipeline_completion/plot_bigwig_files.R")
# DEPENDENCIES: GenomicRanges, rtracklayer
# OUTPUT: svg plots with comparisons for cpm/rpkm/raw and for shifted/raw/deduped for cpm counts
# AUTHOR: LEMR
# DATE: 2025-02-25
# DATE_V1_COMPLETE: 2025-04-08
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
    GENOME_FEATURES <- GenomeInfoDb::keepSeqlevels(GENOME_FEATURES, CHROMOSOME_ROMAN, pruning.mode = "coarse")
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
# Load and process bigwig file
################################################################################
# Setup
EXTENSION_TO_LOAD <- "BW"
EXTENSIONS_TO_REMOVE <- paste0(".", tolower(EXTENSION_TO_LOAD))
ESCAPED_EXTENSIONS <- gsub(".", "\\.", EXTENSIONS_TO_REMOVE, fixed = TRUE)
BIGWIG_FILES <- file_paths_by_type[[EXTENSION_TO_LOAD]]

if (length(BIGWIG_FILES) == 0) {
    stop("No bigwig files found.")
}

#-------------------------------------------------------------------------------
# EXTRACT METADATA FROM FILENAMES
#-------------------------------------------------------------------------------
# Extract base filenames without path and extensions
file_basenames <- basename(BIGWIG_FILES)
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
ADDITIONAL_NUMBER_OF_METADATA_COLUMNS <- 1
metadata_df$bam_type <- sapply(1:nrow(metadata_df), function(row_idx) {
  # Calculate position of BAM type in the row
  bam_type_col_idx <- ORIGINAL_COLUMN_COUNT - ADDITIONAL_NUMBER_OF_METADATA_COLUMNS - metadata_df$none_count[row_idx]
  if (bam_type_col_idx > 0 && bam_type_col_idx <= ORIGINAL_COLUMN_COUNT) {
    return(as.character(metadata_df[row_idx, bam_type_col_idx]))
  } else {
    return("unknown")  # Fallback if position is invalid
  }
})

# Extract the last non-"none" columns normalization
normalization_method <- vector("list", nrow(metadata_df))
for (row_idx in 1:nrow(metadata_df)) {
  end_idx <- ORIGINAL_COLUMN_COUNT - metadata_df$none_count[row_idx]
  if (end_idx >= 2) {
    normalization_method[row_idx] <- metadata_df[row_idx, end_idx]
  }
}

# COMBINE ALL METADATA INTO FINAL DATA FRAME
# Assuming first two columns are sample type and condition idx
final_metadata <- data.frame(
  sample_type = metadata_df[, 1],
  condition_idx = metadata_df[, 2],
  bam_type = metadata_df$bam_type,
  normalization_method = unlist(normalization_method),
  file_path = BIGWIG_FILES,
  stringsAsFactors = FALSE
)

# Assertions
#stopifnot(
#    "No samples are NA." = !any(is.na(final_metadata$sample)),
#    "Metadata has expected dimensions." =
#        nrow(final_metadata) == NUMBER_OF_FILES &&
#        ncol(final_metadata) == length(COLUMN_NAMES),
#    "No duplicate samples." = !any(duplicated(final_metadata))
#)

#################################################################################
# MAIN
#################################################################################
# Create output directory if it doesn't exist
PLOT_OUTPUT_DIR <- file.path(Sys.getenv("HOME"), "data", "preprocessing_test", "plots")
dir.create(PLOT_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

#---------------------------------------------
# Define the comparisons lists
#---------------------------------------------
list_of_comparisons <- list(
  bam_processing_effect = list(
    fixed_columns = c("sample_type", "condition_idx", "normalization_method"),
    vary_column = "bam_type",
    filter = quote(normalization_method == "cpm")
    # Add basename of plot here?
  ),
  normalization_effect = list(
    fixed_columns = c("sample_type", "bam_type", "condition_idx"),
    vary_column = "normalization_method",
    filter = NULL
  ),
  sample_type_comparison = list(
    fixed_columns = c("bam_type", "normalization_method"),
    vary_column = c("sample_type", "condition_idx"),
    filter = NULL
  )
)

message("====================")
for (comparison_name in names(list_of_comparisons)) {
  message("--- Comparison for loop ---")
  comparison <- list_of_comparisons[[comparison_name]]
  fixed_columns_array <- comparison$fixed_columns
  vary_column_value <- comparison$vary_column
  filter_expression <- comparison$filter
  stopifnot(
    "Fixed columns are not in metadata frame. Please adjust." = all(fixed_columns_array %in% names(final_metadata)),
    "Vary column not in metadata frame. Please adjust." = all(vary_column_value %in% names(final_metadata))
  )


  # Apply filter if it exists
  filtered_metadata <- final_metadata
  if (!is.null(filter_expression)) {
    filtered_metadata <- subset(final_metadata, eval(filter_expression, envir = final_metadata))
  }

  message(sprintf("  Dimensions after filtering:%s, %s", nrow(filtered_metadata), ncol(filtered_metadata)))
  # If there are no rows after filtering, skip to next comparison
  if (nrow(filtered_metadata) == 0) {
    message("No data after filtering, skipping comparison")
    next
  }

  # Create a dataframe with just the fixed columns
  fixed_cols_df <- filtered_metadata[, fixed_columns_array, drop = FALSE]

  # Get unique combinations
  unique_combinations <- unique(fixed_cols_df)
  # Process each unique combination
  for (i in 1:nrow(unique_combinations)) {
    message("  --- Unique combination ---")
    # Create a logical vector to match this combination
    row_match <- rep(TRUE, nrow(filtered_metadata))

    for (col in fixed_columns_array) {
      row_match <- row_match & (filtered_metadata[[col]] == unique_combinations[i, col])
    }

    # Subset the dataframe for this combination
    subset_df <- filtered_metadata[row_match, ]

    # Skip if empty subset
    if (nrow(subset_df) == 0) next
    # Create a title describing this comparison
    title_parts <- paste(
      paste(fixed_columns_array, "=",
            sapply(fixed_columns_array, function(col) unique_combinations[i, col])),
      collapse = ", "
    )

    track_title <- paste0(comparison_name, ": ", title_parts)
    print(paste0("  Found ", nrow(subset_df), " tracks for ", track_title))

    # Here you would run your plotting code using subset_df
    # The vary_column(s) will naturally have different values in this subset
    # Optional: Sort the subset by the varying column(s) for consistent track order
    if (length(vary_column_value) == 1) {
      subset_df <- subset_df[order(subset_df[[vary_column_value]]), ]
    }

    # Create a list of fixed values for this combination
    fixed_values <- as.list(unique_combinations[i, ])
    # Start with the comparison name
    filename_parts <- comparison_name
    # Add fixed column values if they exist
    if (length(fixed_values) > 0) {
      # Create sanitized strings for each fixed column
      fixed_strings <- sapply(names(fixed_values), function(col) {
        # Replace spaces and special characters
        value <- gsub("[^a-zA-Z0-9]", "_", fixed_values[[col]])
        return(paste0(col, "-", value))
      })
      # Add to filename parts
      filename_parts <- c(filename_parts, paste(fixed_strings, collapse="_"))
    }
    # Join parts with underscores and add extension
    filename <- paste0(paste(filename_parts, collapse="_"), ".svg")
    # Replace multiple underscores with single one
    filename <- gsub("_+", "_", filename)
    # Create full path
    full_path <- file.path(PLOT_OUTPUT_DIR, filename)
    message(paste0("  File will be saved to ", full_path))

    # [2] Track Container Initialization
    track_container <- vector("list", nrow(subset_df) + 1 + as.numeric(exists("GENOME_FEATURES")))
    track_container[[1]] <- Gviz::GenomeAxisTrack(
        name = sprintf("Chr %s Axis", CHROMOSOME_TO_PLOT),
        fontcolor.title = "black",
        cex.title = 0.7,
        background.title = "white"
    )

    for (file_idx in seq_len(nrow(subset_df))) {
      message("     --- Metadata subsetting ---")
      current_row <- subset_df[file_idx, ]
      current_filepath <- subset_df$file_path[file_idx]

      # Always include sample_type and condition_idx as base information
      base_info <- paste0(
        current_row[["sample_type"]], "_",
        current_row[["condition_idx"]]
      )
      # Add varying column information
      if (length(comparison$vary_column) == 1 &&
          !(comparison$vary_column %in% c("sample_type", "condition_idx"))) {
        # Add single vary column (if it's not already included in base_info)
        vary_col <- comparison$vary_column
        track_name <- paste0(base_info, " (", vary_col, ": ", current_row[[vary_col]], ")")
      } else if (length(comparison$vary_column) > 1) {
        # For multiple vary columns, add any that aren't already in base_info
        additional_cols <- setdiff(comparison$vary_column, c("sample_type", "condition_idx"))
        if (length(additional_cols) > 0) {
          vary_values <- sapply(additional_cols, function(col) {
            paste0(col, "=", current_row[[col]])
          })
          track_name <- paste0(base_info, " (", paste(vary_values, collapse=", "), ")")
        } else {
          track_name <- base_info  # Only sample and condition are varying
        }
      } else {
        # The vary column is already included in base_info
        track_name <- base_info
      }

      message(sprintf("   Importing bigwig file: %s", current_filepath))
      message(paste0("    Adding track: ", track_name))

      bigwig_data <- rtracklayer::import(
        current_filepath,
        format = "BigWig",
        which = GENOME_RANGE_TO_LOAD
      )

      track_container[[file_idx + 1]] <- Gviz::DataTrack(
          range = bigwig_data,
          name = track_name,
          # Apply styling
          showAxis = TRUE,
          showTitle = TRUE,
          type = "h",
          size = 1.2,
          background.title = "white",
          fontcolor.title = "black",
          col.border.title = "#e0e0e0",
          cex.title = 0.7,
          fontface = 1,
          title.width = 1.0,
          col = "darkblue",
          fill = "darkblue"
      )

    }

    message("   All samples imported and added to track_container.")
    # [4] Annotation Track (Conditional)
    if (exists("GENOME_FEATURES")) {
        track_container[[length(track_container)]] <- Gviz::AnnotationTrack(
            GENOME_FEATURES,
            name = "Features",
            size = 0.5,
            background.title = "lightgray",
            fontcolor.title = "black",
            showAxis = FALSE,
            background.panel = "#f5f5f5",
            cex.title = 0.6,
            fill = "#8b4513",
            col = "#8b4513"
        )
    }
    # Create descriptive main plot title
    main_title <- paste0(
      comparison_name, ": ",
      paste(names(fixed_values), "=", unlist(fixed_values), collapse=", ")
    )
    message(sprint("Plotting %s", main_title))

    #svglite::svglite(
    #    filename = full_path,
    #    width = 10,
    #    height = 8,
    #    bg = "white"
    #)

    #Gviz::plotTracks(
    #    trackList = track_container,
    #    chromosome = CHROMOSOME_ROMAN,
    #    from = GENOME_RANGE_TO_LOAD@ranges@start,
    #    to = GENOME_RANGE_TO_LOAD@ranges@width,
    #    margin = 15,
    #    innerMargin = 5,
    #    spacing = 10,
    #    main = main_title,
    #    col.axis = "black",
    #    cex.axis = 0.8,
    #    cex.main = 0.9,
    #    fontface.main = 2,
    #    background.panel = "transparent"
    #)
    #dev.off()
    message("   Plot saved...")

    # Add line break between comparisons for readability
    message("\n")
  }
} # end for loop of comparisons
message("====================")
stop("Breakpoint...")
#---------------------------------------------
# Plot the bigwig processing comparisons
#---------------------------------------------
# Define categories
#UNIQUE_METADATA_CATEGORIES <- list(
#    samples = unique(final_metadata[,"sample"]),
#    bam_processing = unique(final_metadata[,"bam_processing"]),
#    bigwig_processing = unique(final_metadata[, "bigwig_processing"])
#)
#
#VALID_PROCESSING_COMBINATIONS <- expand.grid(
#    UNIQUE_METADATA_CATEGORIES[c("samples", "bam_processing")]
#)
#
## Filter the sample combinations to exclude input shifted
#IS_INPUT_SAMPLE <- VALID_PROCESSING_COMBINATIONS$samples == "input"
#IS_SHIFTED_PROCESSING <- VALID_PROCESSING_COMBINATIONS$bam_processing == "shifted"
#VALID_PROCESSING_COMBINATIONS <- VALID_PROCESSING_COMBINATIONS[!(IS_INPUT_SAMPLE & IS_SHIFTED_PROCESSING), ]

#---------------------------------------------
# Plot the bam processing comparisons for the bigwig files
#---------------------------------------------
# Handle the sample and bam_processing combinations
# Non-printing character to avoid collision
METADATA_COLUMN_SEPARATOR <-  "\x01"

# Create the keys to perform subsetting
columns_for_metadata_keys <- c("sample_type", "condition_idx", "normalization_method")

METADATA_CHARACTER_VECTORS <- lapply(final_metadata[c("sample_type", "bam_processing")], as.character)
METADATA_JOINED_KEYS <- do.call(paste, c(METADATA_CHARACTER_VECTORS, sep = METADATA_COLUMN_SEPARATOR))

#CONTROL_COMBINATIONS_CHARACTERS <- lapply(VALID_PROCESSING_COMBINATIONS, as.character)
#control_joined_keys <- do.call(paste, c(CONTROL_COMBINATIONS_CHARACTERS, sep = METADATA_COLUMN_SEPARATOR))


for (CURRENT_KEY_IDX in seq_along(control_joined_keys)) {
    CURRENT_CONTROL_KEY <- control_joined_keys[CURRENT_KEY_IDX]
    message(sprintf("Processing key %d: %s", CURRENT_KEY_IDX, CURRENT_CONTROL_KEY))

    # [1] Metadata Subsetting
    CURRENT_METADATA_SUBSET <- final_metadata[METADATA_JOINED_KEYS %in% CURRENT_CONTROL_KEY, ]
    message("Metadata subset complete...")
    message("Subset content:")
    print(head(CURRENT_METADATA_SUBSET))

    # [2] Track Container Initialization
    track_container <- vector("list", nrow(CURRENT_METADATA_SUBSET) + 1 + exists("GENOME_FEATURES"))
    track_container[[1]] <- Gviz::GenomeAxisTrack(
        name = sprintf("Chr %s Axis", CHROMOSOME_TO_PLOT),
        fontcolor.title = "black",
        cex.title = 0.7,
        background.title = "white"
    )

    # Process the metadata subet ---------------------
    # [3] Data Track Population
    for (track_idx in seq_len(nrow(CURRENT_METADATA_SUBSET))) {
        track_name <- do.call(paste, c(CURRENT_METADATA_SUBSET[track_idx, 1:3], sep = "_"))
        current_row_filepath <- CURRENT_METADATA_SUBSET[track_idx, "file_paths"]

        message(sprintf("Importing bigwig file: %s", current_row_filepath))
        bigwig_data <- rtracklayer::import(
            current_row_filepath,
            format = "BigWig",
            which = GENOME_RANGE_TO_LOAD
        )

        track_container[[track_idx + 1]] <- Gviz::DataTrack(
            range = bigwig_data,
            name = track_name,
            # Apply styling
            showAxis = TRUE,
            showTitle = TRUE,
            type = "h",
            size = 1.2,
            background.title = "white",
            fontcolor.title = "black",
            col.border.title = "#e0e0e0",
            cex.title = 0.7,
            fontface = 1,
            title.width = 1.0,
            col = "darkblue",
            fill = "darkblue"
        )
        message(sprintf("Sample %s imported...", track_idx))
    } # end for loop of subset metadata ########

    message("All samples imported and added to TRACK_CONTAINER.")

    # [4] Annotation Track (Conditional)
    if (exists("GENOME_FEATURES")) {
        TRACK_CONTAINER[[length(TRACK_CONTAINER)]] <- Gviz::AnnotationTrack(
            GENOME_FEATURES,
            name = "Features",
            size = 0.5,
            background.title = "lightgray",
            fontcolor.title = "black",
            showAxis = FALSE,
            background.panel = "#f5f5f5",
            cex.title = 0.6,
            fill = "#8b4513",
            col = "#8b4513"
        )
    }

    # [5] Plot Generation & Export -----------------
    PLOT_BASENAME <- paste(
        "/bigwig_processing",
        CHROMOSOME_ROMAN,
        gsub(METADATA_COLUMN_SEPARATOR, "_", CURRENT_CONTROL_KEY, fixed = TRUE),
        sep = "_"
    )

    OUTPUT_FILENAME <- paste0(PLOT_OUTPUT_DIR, PLOT_BASENAME, ".svg")
    message(sprintf("Saving file: %s", OUTPUT_FILENAME))

    svglite::svglite(
        filename = OUTPUT_FILENAME,
        width = 10,
        height = 8,
        bg = "white"
    )

    Gviz::plotTracks(
        trackList = TRACK_CONTAINER,
        chromosome = CHROMOSOME_ROMAN,
        from = GENOME_RANGE_TO_LOAD@ranges@start,
        to = GENOME_RANGE_TO_LOAD@ranges@width,
        margin = 15,
        innerMargin = 5,
        spacing = 10,
        col.axis = "black",
        cex.axis = 0.8,
        cex.main = 0.9,
        fontface.main = 2,
        background.panel = "transparent"
    )
    dev.off()
    message(sprintf("Plot saved. Finished processing %d: %s", CURRENT_KEY_IDX, CURRENT_CONTROL_KEY))
} # end keys for loop
message("Finished bigwig processing for loop...")

#---------------------------------------------
# Plot the bam processing comparisons for cpm
#---------------------------------------------
message("Starting processing for bam processing cpm plots...")
# Get the combinations for cpm counts, no need to filter
VALID_BAM_PROCESSING_COMBINATIONS <- expand.grid(
    c(UNIQUE_METADATA_CATEGORIES["samples"],
        bigwig_processing = "cpm")
)

# Create the metadata characters and keys for sample and bigwig comparisons (_SB)
METADATA_CHARACTER_VECTORS_SB <- lapply(final_metadata[c("sample", "bigwig_processing")], as.character)
METADATA_JOINED_KEYS_SB <- do.call(paste, c(METADATA_CHARACTER_VECTORS_SB, sep = METADATA_COLUMN_SEPARATOR))

CONTROL_BAM_COMBINATIONS_CHARACTERS <- lapply(VALID_BAM_PROCESSING_COMBINATIONS, as.character)
control_joined_keys <- do.call(paste, c(CONTROL_BAM_COMBINATIONS_CHARACTERS, sep = METADATA_COLUMN_SEPARATOR))
message("Keys for metadata subsetting created. Starting for loop...")

for (CURRENT_KEY_IDX in seq_along(control_joined_keys)) {
    CURRENT_CONTROL_KEY <- control_joined_keys[CURRENT_KEY_IDX]
    message(sprintf("Processing key %d: %s", CURRENT_KEY_IDX, CURRENT_CONTROL_KEY))

    # [1] Metadata Subsetting
    CURRENT_METADATA_SUBSET <- final_metadata[METADATA_JOINED_KEYS_SB %in% CURRENT_CONTROL_KEY, ]
    message("Metadata subset complete...")
    message("Subset content:")
    print(head(CURRENT_METADATA_SUBSET))

    # [2] Track Container Initialization
    TRACK_CONTAINER <- vector("list", nrow(CURRENT_METADATA_SUBSET) + 1 + exists("GENOME_FEATURES"))
    TRACK_CONTAINER[[1]] <- Gviz::GenomeAxisTrack(
        name = sprintf("Chr %s Axis", CHROMOSOME_TO_PLOT),
        fontcolor.title = "black",
        cex.title = 0.7,
        background.title = "white"
    )

    # Process the metadata subet ---------------------
    # [3] Data Track Population
    for (track_idx in seq_len(nrow(CURRENT_METADATA_SUBSET))) {
        track_name <- do.call(paste, c(CURRENT_METADATA_SUBSET[track_idx, 1:3], sep = "_"))
        current_row_filepath <- CURRENT_METADATA_SUBSET[track_idx, "file_paths"]

        message(sprintf("Importing bigwig file: %s", current_row_filepath))
        bigwig_data <- rtracklayer::import(
            current_row_filepath,
            format = "BigWig",
            which = GENOME_RANGE_TO_LOAD
        )

        TRACK_CONTAINER[[track_idx + 1]] <- Gviz::DataTrack(
            range = bigwig_data,
            name = track_name,
            # Apply styling
            showAxis = TRUE,
            showTitle = TRUE,
            type = "h",
            size = 1.2,
            background.title = "white",
            fontcolor.title = "black",
            col.border.title = "#e0e0e0",
            cex.title = 0.7,
            fontface = 1,
            title.width = 1.0,
            col = "darkblue",
            fill = "darkblue"
        )
        message(sprintf("Sample %s imported...", track_idx))
    } # end for

    message("All samples imported and added to TRACK_CONTAINER.")

    # [4] Annotation Track (Conditional)
    if (exists("GENOME_FEATURES")) {
        TRACK_CONTAINER[[length(TRACK_CONTAINER)]] <- Gviz::AnnotationTrack(
            GENOME_FEATURES,
            name = "Features",
            size = 0.5,
            background.title = "lightgray",
            fontcolor.title = "black",
            showAxis = FALSE,
            background.panel = "#f5f5f5",
            cex.title = 0.6,
            fill = "#8b4513",
            col = "#8b4513"
        )
    }

    # [5] Plot Generation & Export -----------------
    PLOT_BASENAME <- paste(
        "/bigwig_processing",
        CHROMOSOME_ROMAN,
        gsub(METADATA_COLUMN_SEPARATOR, ".", CURRENT_CONTROL_KEY, fixed = TRUE),
        sep = "_"
    )

    OUTPUT_FILENAME <- paste0(PLOT_OUTPUT_DIR, PLOT_BASENAME, ".svg")
    message(sprintf("Saving file: %s", OUTPUT_FILENAME))

    svglite::svglite(
        filename = OUTPUT_FILENAME,
        width = 10,
        height = 8,
        bg = "white"
    )

    Gviz::plotTracks(
        trackList = TRACK_CONTAINER,
        chromosome = CHROMOSOME_ROMAN,
        from = GENOME_RANGE_TO_LOAD@ranges@start,
        to = GENOME_RANGE_TO_LOAD@ranges@width,
        margin = 15,
        innerMargin = 5,
        spacing = 10,
        col.axis = "black",
        cex.axis = 0.8,
        cex.main = 0.9,
        fontface.main = 2,
        background.panel = "transparent"
    )
    dev.off()
    message(sprintf("Plot saved. Finished processing %d: %s", CURRENT_KEY_IDX, CURRENT_CONTROL_KEY))
} # end for loop
message("Finished bam processing for loop...")
