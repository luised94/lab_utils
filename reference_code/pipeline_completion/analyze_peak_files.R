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
CONTROL_SAMPLES <- list(
    "test" = "$HOME/data/241010Bel/alignment/processed_245018_sequence_to_S288C_sorted.bam",
    "input" = "$HOME/data/241010Bel/alignment/processed_245003_sequence_to_S288C_sorted.bam",
    "reference" = "$HOME/data/100303Bel/alignment/processed_034475_sequence_to_S288C_sorted.bam"
)

DATA_DIRECTORY <- file.path(Sys.getenv("HOME"), "data", "preprocessing_test")
SUPPORTED_FILE_TYPES <- c(
  FASTQ = "fastq",
  BAM = "bam",
  BIGWIG = "bw",
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
PEAK_EXTENSIONS <- c("NARROW_PEAK", "BED", "XLS", "GAPPED_PEAK", "BROAD_PEAK")
PEAK_FILES <- file_paths_by_type[[PEAK_EXTENSIONS]]
NUMBER_OF_FILES <- length(PEAK_FILES)
# Load PEAK_FILE_COLUMNS variable
source(normalizePath("~/lab_utils/core_scripts/peak_file_columns_by_type.R"))
stopifnot("Problem loading PEAK_FILE_COLUMNS variable from R file." = exists(PEAK_FILE_COLUMNS))
FILE_EXTENSION_TO_REMOVE <- paste0(".", SUPPORTED_FILE_TYPES)
EXPECTED_NUM_UNDERSCORES <- 6
#
file_basenames <- basename(BIGWIG_FILES)
filenames_without_extension <- unlist(
  lapply(file_basenames, function(file_basename, extensions_to_remove){
    for (EXTENSION_IDX in 1:length(extensions_to_remove)) {
      extension_to_remove <- extensions_to_remove[EXTENSION_IDX]
      if (grepl(pattern = extension_to_remove, x = file_basename)) {
        filename_without_extension <- gsub(pattern = extension_to_remove, replacement = "", x = file_basename)
      }
    }
    return(filename_without_extension)
  }, extensions_to_remove = FILE_EXTENSION_TO_REMOVE)
)
split_metadata <- strsplit(filenames_without_extension, split = "_")

underscore_counts <- lengths(gregexpr("_", filenames_without_extension))
NUMBER_OF_COLUMNS <- lengths(PEAK_FILE_COLUMNS)

if (length(unique(underscore_counts)) != 1 || unique(underscore_counts) != EXPECTED_NUM_UNDERSCORES) {
  stop("Filenames must contain exactly ", EXPECTED_NUM_UNDERSCORES, " underscores.")
}

## Preallocation
#bigwig_metadata_df <- data.frame(matrix(NA, nrow = NUMBER_OF_FILES, ncol = NUMBER_OF_COLUMNS))
#
#for (row_index in 1:nrow(bigwig_metadata_df)) {
#    bigwig_metadata_df[row_index, 1:3] <- split_metadata[[row_index]]
#}
#
#colnames(bigwig_metadata_df) <- COLUMN_NAMES
#bigwig_metadata_df$file_paths <- BIGWIG_FILES
#
## Assertions
#stopifnot(
#    "No samples are NA." = !any(is.na(bigwig_metadata_df$sample)),
#    "Metadata has expected dimensions." =
#        nrow(bigwig_metadata_df) == NUMBER_OF_FILES &&
#        ncol(bigwig_metadata_df) == length(COLUMN_NAMES),
#    "No duplicate samples." = !any(duplicated(bigwig_metadata_df))
#)
#
#if (length(BIGWIG_FILES) == 0) {
#    stop("No bigwig files found.")
#}
#
##################################################################################
## MAIN
##################################################################################
##---------------------------------------------
## Plot the bigwig processing comparisons
##---------------------------------------------
## Define categories
#UNIQUE_METADATA_CATEGORIES <- list(
#    samples = unique(bigwig_metadata_df[,"sample"]),
#    bam_processing = unique(bigwig_metadata_df[,"bam_processing"]),
#    bigwig_processing = unique(bigwig_metadata_df[, "bigwig_processing"])
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
#
## Handle the sample and bam_processing combinations
## Non-printing character to avoid collision
#METADATA_COLUMN_SEPARATOR <-  "\x01"
#
## Create the keys to perform subsetting
#METADATA_CHARACTER_VECTORS <- lapply(bigwig_metadata_df[c("sample", "bam_processing")], as.character)
#METADATA_JOINED_KEYS <- do.call(paste, c(METADATA_CHARACTER_VECTORS, sep = METADATA_COLUMN_SEPARATOR))
#
#CONTROL_COMBINATIONS_CHARACTERS <- lapply(VALID_PROCESSING_COMBINATIONS, as.character)
#control_joined_keys <- do.call(paste, c(CONTROL_COMBINATIONS_CHARACTERS, sep = METADATA_COLUMN_SEPARATOR))
#
## Create output directory if it doesn't exist
#PLOT_OUTPUT_DIR <- file.path(Sys.getenv("HOME"), "data", "preprocessing_test", "plots")
#dir.create(PLOT_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
#
#for (CURRENT_KEY_IDX in seq_along(control_joined_keys)) {
#    CURRENT_CONTROL_KEY <- control_joined_keys[CURRENT_KEY_IDX]
#    message(sprintf("Processing key %d: %s", CURRENT_KEY_IDX, CURRENT_CONTROL_KEY))
#
#    # [1] Metadata Subsetting
#    CURRENT_METADATA_SUBSET <- bigwig_metadata_df[METADATA_JOINED_KEYS %in% CURRENT_CONTROL_KEY, ]
#    message("Metadata subset complete...")
#    message("Subset content:")
#    print(head(CURRENT_METADATA_SUBSET))
#
#    # [2] Track Container Initialization
#    TRACK_CONTAINER <- vector("list", nrow(CURRENT_METADATA_SUBSET) + 1 + exists("GENOME_FEATURES"))
#    TRACK_CONTAINER[[1]] <- Gviz::GenomeAxisTrack(
#        name = sprintf("Chr %s Axis", CHROMOSOME_TO_PLOT),
#        fontcolor.title = "black",
#        cex.title = 0.7,
#        background.title = "white"
#    )
#
#    # Process the metadata subet ---------------------
#    # [3] Data Track Population
#    for (track_idx in seq_len(nrow(CURRENT_METADATA_SUBSET))) {
#        track_name <- do.call(paste, c(CURRENT_METADATA_SUBSET[track_idx, 1:3], sep = "_"))
#        current_row_filepath <- CURRENT_METADATA_SUBSET[track_idx, "file_paths"]
#
#        message(sprintf("Importing bigwig file: %s", current_row_filepath))
#        bigwig_data <- rtracklayer::import(
#            current_row_filepath,
#            format = "BigWig",
#            which = GENOME_RANGE_TO_LOAD
#        )
#
#        TRACK_CONTAINER[[track_idx + 1]] <- Gviz::DataTrack(
#            range = bigwig_data,
#            name = track_name,
#            # Apply styling
#            showAxis = TRUE,
#            showTitle = TRUE,
#            type = "h",
#            size = 1.2,
#            background.title = "white",
#            fontcolor.title = "black",
#            col.border.title = "#e0e0e0",
#            cex.title = 0.7,
#            fontface = 1,
#            title.width = 1.0,
#            col = "darkblue",
#            fill = "darkblue"
#        )
#        message(sprintf("Sample %s imported...", track_idx))
#    } # end for
#
#    message("All samples imported and added to TRACK_CONTAINER.")
#
#    # [4] Annotation Track (Conditional)
#    if (exists("GENOME_FEATURES")) {
#        TRACK_CONTAINER[[length(TRACK_CONTAINER)]] <- Gviz::AnnotationTrack(
#            GENOME_FEATURES,
#            name = "Features",
#            size = 0.5,
#            background.title = "lightgray",
#            fontcolor.title = "black",
#            showAxis = FALSE,
#            background.panel = "#f5f5f5",
#            cex.title = 0.6,
#            fill = "#8b4513",
#            col = "#8b4513"
#        )
#    }
#
#    # [5] Plot Generation & Export -----------------
#    PLOT_BASENAME <- paste(
#        "/bigwig_processing",
#        CHROMOSOME_ROMAN,
#        gsub(METADATA_COLUMN_SEPARATOR, "_", CURRENT_CONTROL_KEY, fixed = TRUE),
#        sep = "_"
#    )
#
#    OUTPUT_FILENAME <- paste0(PLOT_OUTPUT_DIR, PLOT_BASENAME, ".svg")
#    message(sprintf("Saving file: %s", OUTPUT_FILENAME))
#
#    svglite::svglite(
#        filename = OUTPUT_FILENAME,
#        width = 10,
#        height = 8,
#        bg = "white"
#    )
#
#    Gviz::plotTracks(
#        trackList = TRACK_CONTAINER,
#        chromosome = CHROMOSOME_ROMAN,
#        from = GENOME_RANGE_TO_LOAD@ranges@start,
#        to = GENOME_RANGE_TO_LOAD@ranges@width,
#        margin = 15,
#        innerMargin = 5,
#        spacing = 10,
#        col.axis = "black",
#        cex.axis = 0.8,
#        cex.main = 0.9,
#        fontface.main = 2,
#        background.panel = "transparent"
#    )
#    dev.off()
#    message(sprintf("Plot saved. Finished processing %d: %s", CURRENT_KEY_IDX, CURRENT_CONTROL_KEY))
#} # end for loop
#message("Finished bigwig processing for loop...")
#
##---------------------------------------------
## Plot the bam processing comparisons for cpm
##---------------------------------------------
#message("Starting processing for bam processing cpm plots...")
## Get the combinations for cpm counts, no need to filter
#VALID_BAM_PROCESSING_COMBINATIONS <- expand.grid(
#    c(UNIQUE_METADATA_CATEGORIES["samples"],
#        bigwig_processing = "cpm")
#)
#
## Create the metadata characters and keys for sample and bigwig comparisons (_SB)
#METADATA_CHARACTER_VECTORS_SB <- lapply(bigwig_metadata_df[c("sample", "bigwig_processing")], as.character)
#METADATA_JOINED_KEYS_SB <- do.call(paste, c(METADATA_CHARACTER_VECTORS_SB, sep = METADATA_COLUMN_SEPARATOR))
#
#CONTROL_BAM_COMBINATIONS_CHARACTERS <- lapply(VALID_BAM_PROCESSING_COMBINATIONS, as.character)
#control_joined_keys <- do.call(paste, c(CONTROL_BAM_COMBINATIONS_CHARACTERS, sep = METADATA_COLUMN_SEPARATOR))
#message("Keys for metadata subsetting created. Starting for loop...")
#
#for (CURRENT_KEY_IDX in seq_along(control_joined_keys)) {
#    CURRENT_CONTROL_KEY <- control_joined_keys[CURRENT_KEY_IDX]
#    message(sprintf("Processing key %d: %s", CURRENT_KEY_IDX, CURRENT_CONTROL_KEY))
#
#    # [1] Metadata Subsetting
#    CURRENT_METADATA_SUBSET <- bigwig_metadata_df[METADATA_JOINED_KEYS_SB %in% CURRENT_CONTROL_KEY, ]
#    message("Metadata subset complete...")
#    message("Subset content:")
#    print(head(CURRENT_METADATA_SUBSET))
#
#    # [2] Track Container Initialization
#    TRACK_CONTAINER <- vector("list", nrow(CURRENT_METADATA_SUBSET) + 1 + exists("GENOME_FEATURES"))
#    TRACK_CONTAINER[[1]] <- Gviz::GenomeAxisTrack(
#        name = sprintf("Chr %s Axis", CHROMOSOME_TO_PLOT),
#        fontcolor.title = "black",
#        cex.title = 0.7,
#        background.title = "white"
#    )
#
#    # Process the metadata subet ---------------------
#    # [3] Data Track Population
#    for (track_idx in seq_len(nrow(CURRENT_METADATA_SUBSET))) {
#        track_name <- do.call(paste, c(CURRENT_METADATA_SUBSET[track_idx, 1:3], sep = "_"))
#        current_row_filepath <- CURRENT_METADATA_SUBSET[track_idx, "file_paths"]
#
#        message(sprintf("Importing bigwig file: %s", current_row_filepath))
#        bigwig_data <- rtracklayer::import(
#            current_row_filepath,
#            format = "BigWig",
#            which = GENOME_RANGE_TO_LOAD
#        )
#
#        TRACK_CONTAINER[[track_idx + 1]] <- Gviz::DataTrack(
#            range = bigwig_data,
#            name = track_name,
#            # Apply styling
#            showAxis = TRUE,
#            showTitle = TRUE,
#            type = "h",
#            size = 1.2,
#            background.title = "white",
#            fontcolor.title = "black",
#            col.border.title = "#e0e0e0",
#            cex.title = 0.7,
#            fontface = 1,
#            title.width = 1.0,
#            col = "darkblue",
#            fill = "darkblue"
#        )
#        message(sprintf("Sample %s imported...", track_idx))
#    } # end for
#
#    message("All samples imported and added to TRACK_CONTAINER.")
#
#    # [4] Annotation Track (Conditional)
#    if (exists("GENOME_FEATURES")) {
#        TRACK_CONTAINER[[length(TRACK_CONTAINER)]] <- Gviz::AnnotationTrack(
#            GENOME_FEATURES,
#            name = "Features",
#            size = 0.5,
#            background.title = "lightgray",
#            fontcolor.title = "black",
#            showAxis = FALSE,
#            background.panel = "#f5f5f5",
#            cex.title = 0.6,
#            fill = "#8b4513",
#            col = "#8b4513"
#        )
#    }
#
#    # [5] Plot Generation & Export -----------------
#    PLOT_BASENAME <- paste(
#        "/bigwig_processing",
#        CHROMOSOME_ROMAN,
#        gsub(METADATA_COLUMN_SEPARATOR, ".", CURRENT_CONTROL_KEY, fixed = TRUE),
#        sep = "_"
#    )
#
#    OUTPUT_FILENAME <- paste0(PLOT_OUTPUT_DIR, PLOT_BASENAME, ".svg")
#    message(sprintf("Saving file: %s", OUTPUT_FILENAME))
#
#    svglite::svglite(
#        filename = OUTPUT_FILENAME,
#        width = 10,
#        height = 8,
#        bg = "white"
#    )
#
#    Gviz::plotTracks(
#        trackList = TRACK_CONTAINER,
#        chromosome = CHROMOSOME_ROMAN,
#        from = GENOME_RANGE_TO_LOAD@ranges@start,
#        to = GENOME_RANGE_TO_LOAD@ranges@width,
#        margin = 15,
#        innerMargin = 5,
#        spacing = 10,
#        col.axis = "black",
#        cex.axis = 0.8,
#        cex.main = 0.9,
#        fontface.main = 2,
#        background.panel = "transparent"
#    )
#    dev.off()
#    message(sprintf("Plot saved. Finished processing %d: %s", CURRENT_KEY_IDX, CURRENT_CONTROL_KEY))
#} # end for loop
#message("Finished bam processing for loop...")
