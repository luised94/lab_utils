###############################################################################
# Plot bigwig files
################################################################################
# PURPOSE: Plot bigwig files generated during the pipeline completion tests
# Conclusion: {{fill}}
# USAGE: source("reference_code/pipeline_completion/plot_bigwig_files.R")
# DEPENDENCIES: GenomicRanges, rtracklayer
# OUTPUT: {{fill}}
# AUTHOR: LEMR
# DATE: 2025-02-25
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
BIGWIG_FILES <- file_paths_by_type[["BIGWIG"]]
NUMBER_OF_FILES <- length(BIGWIG_FILES)
COLUMN_NAMES <- c("sample", "bam_processing", "bigwig_processing", "file_paths")
FILE_EXTENSION_TO_REMOVE <- ".bw"
EXPECTED_NUM_UNDERSCORES <- 2

file_basenames <- basename(BIGWIG_FILES)
filenames_without_extension <- gsub(FILE_EXTENSION_TO_REMOVE, "", file_basenames)
split_metadata <- strsplit(filenames_without_extension, split = "_")

underscore_counts <- lengths(gregexpr("_", filenames_without_extension))
NUMBER_OF_COLUMNS <- length(COLUMN_NAMES)

if (length(unique(underscore_counts)) != 1 || unique(underscore_counts) != EXPECTED_NUM_UNDERSCORES) {
  stop("Filenames must contain exactly ", EXPECTED_NUM_UNDERSCORES, " underscores.")
}

# Preallocation
bigwig_metadata_df <- data.frame(matrix(NA, nrow = NUMBER_OF_FILES, ncol = NUMBER_OF_COLUMNS))

for (row_index in 1:nrow(bigwig_metadata_df)) {
    bigwig_metadata_df[row_index, 1:3] <- split_metadata[[row_index]]
}

colnames(bigwig_metadata_df) <- COLUMN_NAMES
bigwig_metadata_df$file_paths <- BIGWIG_FILES

# Assertions
stopifnot(
    "No samples are NA." = !any(is.na(bigwig_metadata_df$sample)),
    "Metadata has expected dimensions." =
        nrow(bigwig_metadata_df) == NUMBER_OF_FILES &&
        ncol(bigwig_metadata_df) == length(COLUMN_NAMES),
    "No duplicate samples." = !any(duplicated(bigwig_metadata_df))
)

if (length(BIGWIG_FILES) == 0) {
    stop("No bigwig files found.")
}

#################################################################################
# MAIN
#################################################################################
# Define track style defaults
TRACK_STYLE <- list(
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
    ylim = c(0, 1000),
    col = "darkblue",
    fill = "darkblue"
)

# Define plot style defaults
PLOT_STYLE <- list(
    margin = 15,
    innerMargin = 5,
    spacing = 10,
    col.axis = "black",
    cex.axis = 0.8,
    cex.main = 0.9,
    fontface.main = 2,
    background.panel = "transparent"
)

# Define categories
CATEGORIES <- list(
    samples = unique(bigwig_metadata_df[,"sample"]),
    bam_processing = unique(bigwig_metadata_df[,"bam_processing"]),
    bigwig_processing = unique(bigwig_metadata_df[, "bigwig_processing"])
)

sample_bam_combinations <- expand.grid(
    CATEGORIES[c("samples", "bam_processing")]
)

# Filter the sample combinations to exclude input shifted
is_input_combination <- sample_bam_combinations$samples == "input" 
is_shifted_combination <- sample_bam_combinations$bam_processing == "shifted"
sample_bam_combinations <- sample_bam_combinations[!(is_input_combination & is_shifted_combination), ]


# Handle the sample and bam_processing combinations
COLUMN_SEPARATOR <-  "\x01"
metadata_chr <- lapply(bigwig_metadata_df[c("sample", "bam_processing")], as.character)
metadata_keys <- do.call(paste, c(metadata_chr, sep = COLUMN_SEPARATOR))

control_combos_chr <- lapply(sample_bam_combinations, as.character)
control_keys <- do.call(paste, c(control_combos_chr, sep = COLUMN_SEPARATOR))

for (key_idx in 1:length(control_keys)) {
    control_key <- control_keys[key_idx]
    keys_in_current_key <- metadata_keys %in% control_key
    cat(sprintf("Processing %s key\n", key_idx))
    cat(sprintf("Using key = %s \n", control_key))


    rows_to_analyze <- bigwig_metadata_df[keys_in_current_key, ]
    message("Sample subsetting complete...")
    print(head(rows_to_analyze))

    track_list <- vector("list", nrow(rows_to_analyze) + 1 + 1)
    track_list[[1]] <- Gviz::GenomeAxisTrack(
        name = sprintf("Chr %s Axis", CHROMOSOME_TO_PLOT)
    )

    for (i in 1:nrow(rows_to_analyze)) {
        track_name <- do.call(paste, c(rows_to_analyze[i, 1:3], sep = "_"))
        rows_file_path <- rows_to_analyze[i, "file_paths"]
        message(sprintf("Importing file path %s...", rows_file_path))
        bigwig_data <- rtracklayer::import(
            rows_file_path,
            format = "BigWig",
            which = GENOME_RANGE_TO_LOAD
        )
        track_list[[i + 1]] <- Gviz::DataTrack(
            data = bigwig_data,
            name = track_name
        )
        message(sprintf("Sample %s imported...", i))
    }

    message("All samples imported and added to track_list.")

    if (exists("GENOME_FEATURES")) {
        track_list[[length(track_list)]] <- Gviz::AnnotationTrack(
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

    # Style the tracks
    # Special styling for axis track
    Gviz::displayPars(track_list[[1]]) <- list(
        fontcolor.title = "black",
        cex.title = 0.7,
        background.title = "white"
    )

    # Apply styles to tracks and plot
    # Skip axis track and feature track
    for (i in 2:(length(track_list)-1)) {
        Gviz::displayPars(track_list[[i]]) <- TRACK_STYLE
    }

    Gviz::plotTracks(
        trackList = track_list,
        PLOT_STYLE
    )
}

#
#sample_bigwig_combinations <- expand.grid(
#    c(CATEGORIES["samples"],
#        bigwig_processing = "cpm")
#)
#metadata_chr <- lapply(bigwig_metadata_df[c("sample", "bigwig_processing")], as.character)
#metadata_keys <- do.call(paste, c(metadata_chr, sep = COLUMN_SEPARATOR))
#
#control_combos_chr <- lapply(sample_bigwig_combinations, as.character)
#control_keys <- do.call(paste, c(control_combos_chr, sep = COLUMN_SEPARATOR))
#is_not_input <- !grepl("input", control_keys)
#control_keys <- control_keys[is_not_input]
# Determine preallocation parameters for tracks list
## Initialize track list with genome axis
#tracks <- list(
#    Gviz::GenomeAxisTrack(
#        name = sprintf("Chr %s Axis", CHROMOSOME_TO_PLOT)
#    )
#)
#
#for (bigwig_file in BIGWIG_FILES) {
#    bigwig_data <- rtracklayer::import(bigwig_file_path, which = genomic_range)
#    genome_track <- Gviz::DataTrack(bigwig_data)
#    tracks[[length(tracks) + 1]] <- genome_track
#}
#if (exists("features")) {
#    tracks[[length(tracks) + 1]] <- Gviz::AnnotationTrack(
#        features,
#        name = "Features",
#        size = 0.5,
#        background.title = "lightgray",
#        fontcolor.title = "black",
#        cex.title = 0.6
#    )
#}
#
#OUTPUT_FILE_NAME
#svglite(
#    filename = OUTPUT_FILE_NAME,
#    width = 10,
#    height = 8,
#    bg = "white",
#    pointsize = 12,
#    standalone = TRUE,
#    system_fonts = list(),
#    user_fonts = list(),
#    web_fonts = list(),
#    id = NULL,
#    fix_text_size = TRUE,
#    scaling = 1,
#    always_valid = FALSE,
#    file
#)
#Gviz::plotTracks(
#    trackList = tracks
#)
#dev.off()
