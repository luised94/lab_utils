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
current_timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
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
final_metadata_df <- data.frame(
  sample_type = metadata_df[, 1],
  condition_idx = metadata_df[, 2],
  bam_type = metadata_df$bam_type,
  normalization_method = unlist(normalization_method),
  file_path = BIGWIG_FILES,
  stringsAsFactors = FALSE
)

# Determine unique values for each categories
message("Determining unique categories...")
metadata_columns <- setdiff(names(final_metadata_df), "file_path")
unique_categories_lst <- vector("list", length(metadata_columns))
message("Before...")
print(unique_categories_lst)
for (array_index in 1:length(metadata_columns)) {
  column_name_chr <- metadata_columns[array_index]
  print(column_name_chr)
  unique_values <- unique(final_metadata_df[, column_name_chr])
  print(unique_values)
  unique_categories_lst[[array_index]] <- unique_values
}
message("After...")
print(unique_categories_lst)
names(unique_categories_lst) <- metadata_columns
lapply(names(unique_categories_lst), function(column_name_chr){
  message(sprintf("--- Column %s ---", column_name_chr))
  #message("  ---Style 1 ---")
  #for (value in unique_categories_lst[[column_name_chr]]) {
    #pad_value <- paste0("|_ ", value)
    #print(pad_value)
  #}
  message("  --- Style 2 ---")
  collapsed_values <- paste(unique_categories_lst[[column_name_chr]], collapse = ",")
  print(collapsed_values)
  return()
})
stop("Breakpoint...")

# Assertions
#stopifnot(
#    "No samples are NA." = !any(is.na(final_metadata_df$sample)),
#    "Metadata has expected dimensions." =
#        nrow(final_metadata_df) == NUMBER_OF_FILES &&
#        ncol(final_metadata_df) == length(COLUMN_NAMES),
#    "No duplicate samples." = !any(duplicated(final_metadata_df))
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
title_comparison_template = paste(
  "Comparison: %s",   # Comparison ID
  "Chromosome %s", # Chr info
  "Date: %s",
  "Vary Col(s): %s",
  "Fixed Col(s): %s",
  "Fixed values: %s",
  sep = "\n"
)
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
# TODO: Add scaling mode
# scaling_modes <- c("local", "individual")
message("====================")
for (comparison_name in names(list_of_comparisons)) {
  message("--- Comparison for loop ---")
  comparison <- list_of_comparisons[[comparison_name]]
  fixed_columns_array <- comparison$fixed_columns
  vary_column_value <- comparison$vary_column
  filter_expression <- comparison$filter
  stopifnot(
    "Fixed columns are not in metadata frame. Please adjust." = all(fixed_columns_array %in% names(final_metadata_df)),
    "Vary column not in metadata frame. Please adjust." = all(vary_column_value %in% names(final_metadata_df))
  )

  # Apply filter if it exists
  filtered_metadata_df <- final_metadata_df
  if (!is.null(filter_expression)) {
    filtered_metadata_df <- subset(final_metadata_df, eval(filter_expression, envir = final_metadata_df))
  }

  message(sprintf("  Dimensions after filtering:%s, %s", nrow(filtered_metadata_df), ncol(filtered_metadata_df)))
  # If there are no rows after filtering, skip to next comparison
  if (nrow(filtered_metadata_df) == 0) {
    message("No data after filtering, skipping comparison")
    next
  }

  # Create a dataframe with just the fixed columns
  fixed_cols_df <- filtered_metadata_df[, fixed_columns_array, drop = FALSE]

  # Get unique combinations
  unique_combinations_df <- unique(fixed_cols_df)
  # Process each unique combination
  for (i in 1:nrow(unique_combinations_df)[1]) {
    message("  --- Unique combination ---")
    message(sprintf("Processing idx %s ------", i))
    # Create a logical vector to match this combination
    row_match <- rep(TRUE, nrow(filtered_metadata_df))

    for (col in fixed_columns_array) {
      row_match <- row_match & (filtered_metadata_df[[col]] == unique_combinations_df[i, col])
    }

    # Subset the dataframe for this combination
    subset_df <- filtered_metadata_df[row_match, ]

    # Skip if empty subset
    if (nrow(subset_df) == 0) next
    # Create a title describing this comparison
    comparison_description_chr <- paste(
      paste(fixed_columns_array, "=",
            sapply(fixed_columns_array, function(col) unique_combinations_df[i, col])),
      collapse = ", "
    )
    print(paste0("  Found ", nrow(subset_df), " tracks for ", comparison_description_chr))

    # Here you would run your plotting code using subset_df
    # The vary_column(s) will naturally have different values in this subset
    # Optional: Sort the subset by the varying column(s) for consistent track order
    if (length(vary_column_value) == 1) {
      subset_df <- subset_df[order(subset_df[[vary_column_value]]), ]
    }

    # Create a list of fixed values for this combination
    fixed_values_list <- as.list(unique_combinations_df[i, ])
    # Start with the comparison name
    filename_parts <- comparison_name
    # Add fixed column values if they exist
    if (length(fixed_values_list) > 0) {
      # Create sanitized strings for each fixed column
      fixed_strings <- sapply(names(fixed_values_list), function(col) {
        # Replace spaces and special characters
        value <- gsub("[^a-zA-Z0-9]", "_", fixed_values_list[[col]])
        return(paste0(col, "-", value))
      })
      # Add to filename parts
      filename_parts <- c(filename_parts, paste(fixed_strings, collapse="_"))
    }
    # Join parts with underscores and add extension
    filename <- paste0(paste(filename_parts, collapse="_"), ".svg")
    # Replace multiple underscores with single one
    filename <- gsub("_+", "_", filename)
    filename <- paste0(current_timestamp, "_", filename)

    # Create full path
    plot_output_path <- file.path(PLOT_OUTPUT_DIR, filename)
    message(paste0("  File will be saved to ", plot_output_path))

    # TODO: Will break with current_timestamp prefix appending //
    # current_timestamp is always unique essentially. //
    # Determine files in directory up front. Parse the plot output path and match. //
    if(file.exists(plot_output_path)) {
      message(sprintf("[DUPLICATE] File %s already exists. Skipping..."))
      next
    }

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
        track_name <- paste0(base_info, "_", current_row[[vary_col]])
        #track_name <- paste0(base_info, " (", vary_col, ": ", current_row[[vary_col]], ")")
      } else if (length(comparison$vary_column) > 1) {
        # For multiple vary columns, add any that aren't already in base_info
        additional_cols <- setdiff(comparison$vary_column, c("sample_type", "condition_idx"))
        if (length(additional_cols) > 0) {
          vary_values <- sapply(additional_cols, function(col) {
            current_row[[col]]
            #paste0(col, ":", current_row[[col]])
          })
          track_name <- paste0(base_info, paste(vary_values, collapse="_"))
        } else {
          track_name <- base_info  # Only sample and condition are varying
        }
      } else {
        # The vary column is already included in base_info
        track_name <- base_info
      }

      message(sprintf("    Importing bigwig file: %s", current_filepath))
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
    plot_title_chr <- sprintf(
        title_comparison_template,
        comparison_name,
        CHROMOSOME_TO_PLOT,
        current_timestamp,
        paste(vary_column_value, collapse = ","),
        paste(fixed_columns_array, collapse = ","),
        paste(unlist(fixed_values_list), collapse = ",")
    )
    message("~~~~~~~~~~~~~~~~~~~~~~")
    message("Plotting...")
    #print(title_comparison_template)
    #print(comparison_name)
    #print(CHROMOSOME_TO_PLOT)
    #print(current_timestamp)
    #print(paste(vary_column_value, collapse = ","))
    #print(paste(fixed_columns_array, collapse = ","))
    #print(paste(unlist(fixed_values_list), collapse = ","))
    message(plot_title_chr)
    message("~~~~~~~~~~~~~~~~~~~~~~")

    svglite::svglite(
        filename = plot_output_path,
        width = 10,
        height = 8,
        bg = "white"
    )

    Gviz::plotTracks(
        trackList = track_container,
        chromosome = CHROMOSOME_ROMAN,
       from = GENOME_RANGE_TO_LOAD@ranges@start,
        to = GENOME_RANGE_TO_LOAD@ranges@width,
        margin = 15,
        innerMargin = 5,
        spacing = 10,
        main = plot_title_chr,
        col.axis = "black",
        cex.axis = 0.8,
        cex.main = 0.9,
       fontface.main = 2,
        background.panel = "transparent"
    )
    dev.off()
    message("   Plot saved...")
    break
    # Add line break between comparisons for readability
    message("\n")
  } # end for loop for unique combinations of each comparison
  break
} # end for loop of comparisons
message("====================")
message("Finished bigwig processing for loop...")
