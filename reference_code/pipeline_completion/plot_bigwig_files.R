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
FILE_FEATURE_DIRECTORY <- file.path(Sys.getenv("HOME"), "data", "feature_files")
stopifnot(
    "Genome directory not found" = dir.exists(FILE_GENOME_DIRECTORY),
    "Feature directory not found" = dir.exists(FILE_FEATURE_DIRECTORY)
)

# Load reference genome
REF_GENOME_FILE <- list.files(
    path = FILE_GENOME_DIRECTORY,
    pattern = FILE_FEATURE_DIRECTORY,
    full.names = TRUE,
    recursive = TRUE
)[1]

#if (length(REF_GENOME_FILE) == 0) {
#    stop(sprintf("No reference genome files found matching pattern '%s' in: %s",
#                FILE_GENOME_DIRECTORY,
#                FILE_FEATURE_DIRECTORY))
#}
#
#if (!file.exists(REF_GENOME_FILE)) {
#    stop(sprintf("Reference genome file not accessible: %s", REF_GENOME_FILE[1]))
#}

#REFERENCE_GENOME_DSS <- Biostrings::readDNAStringSet(REF_GENOME_FILE)
#
## Create chromosome range
#CHROMOSOME_TO_PLOT <- RUNTIME_CONFIG$process_chromosome
#CHROMOSOME_WIDTH <- REFERENCE_GENOME_DSS[CHROMOSOME_TO_PLOT]@ranges@width
#CHROMOSOME_ROMAN <- paste0("chr", utils::as.roman(CHROMOSOME_TO_PLOT))
#
#GENOME_RANGE_TO_LOAD <- GenomicRanges::GRanges(
#    seqnames = CHROMOSOME_ROMAN,
#    ranges = IRanges::IRanges(start = 1, end = CHROMOSOME_WIDTH),
#    strand = "*"
#)
#
## Load feature file (annotation)
#FEATURE_FILE <- list.files(
#    path = GENOME_TRACK_CONFIG$file_feature_directory,
#    pattern = GENOME_TRACK_CONFIG$file_feature_pattern,
#    full.names = TRUE
#)[1]
#
#if (length(FEATURE_FILE) == 0)
#{
#    warning(sprintf("No feature files found matching pattern '%s' in: %s",
#                   GENOME_TRACK_CONFIG$file_feature_pattern,
#                   GENOME_TRACK_CONFIG$file_feature_directory))
#}
#
#if (!is.null(FEATURE_FILE))
#{
#    GENOME_FEATURES <- rtracklayer::import(FEATURE_FILE)
#    # Convert to chrRoman format
#    GenomeInfoDb::seqlevels(GENOME_FEATURES) <- paste0(
#        "chr",
#        utils::as.roman(gsub("chr", "", GenomeInfoDb::seqlevels(GENOME_FEATURES)))
#    )
#}
#
################################################################################
# Setup directories and file paths
################################################################################
CONTROL_SAMPLES <- list(
    "test" = "$HOME/data/241010Bel/alignment/processed_245018_sequence_to_S288C_sorted.bam",
    "input" = "$HOME/data/241010Bel/alignment/processed_245003_sequence_to_S288C_sorted.bam",
    "reference" = "$HOME/data/100303Bel/alignment/processed_034475_sequence_to_S288C_sorted.bam"
)

DATA_DIRECTORY <- file.path(Sys.getenv("HOME"), "preprocessing_test")
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
    file_pattern <- paste0("\\.", extension, "$")
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


