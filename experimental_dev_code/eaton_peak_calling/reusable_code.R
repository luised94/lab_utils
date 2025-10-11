

######### JUNK ###########
#CHROMOSOME_TO_PLOT <- 10
#FILE_FEATURE_DIRECTORY <- file.path(Sys.getenv("HOME"), "data", "feature_files")
#FILE_GENOME_DIRECTORY <- file.path(Sys.getenv("HOME"), "data", "REFGENS")
#
#FILE_GENOME_PATTERN <- "S288C_refgenome.fna"
#FILE_FEATURE_PATTERN <- "eaton_peaks"
## @ques: maybe I should download this file locally.
#REFERENCE_BAM_PATH <- "/home/luised94/data/100303Bel/alignment/consolidated_034475_sequence_to_S288C_sorted.bam"
#REFERENCE_BAM_HEADER <- Rsamtools::scanBamHeader(REFERENCE_BAM_PATH)
#print(CHROMOSOME_LENGTHS_nls)
#
#stopifnot(
#  #"REFERENCE_BAM_PATH does not exist." =
#  #file.exists(REFERENCE_BAM_PATH),
#  "Genome directory not found" =
#  dir.exists(FILE_GENOME_DIRECTORY),
#  "Feature directory not found" =
#  dir.exists(FILE_FEATURE_DIRECTORY),
#  "CHROMOSOME_TO_PLOT must be numeric." =
#  is.numeric(CHROMOSOME_TO_PLOT)
#)
#
#required_packages <- c(
#  "rtracklayer",
#  "GenomicRanges",
#  "Gviz"
#)
#
#is_missing_package <- !sapply(
#  X = required_packages,
#  FUN = requireNamespace,
#  quietly = TRUE
#)
#
#if (
#  !is.character(required_packages) ||
#  length(required_packages) == 0
#) {
#  stop("required_packages must be a non-empty character vector")
#}
#
#missing_packages <- required_packages[is_missing_package]
#if (length(missing_packages) > 0 ) {
#  stop(
#  "Missing packages. Please install using renv:\n",
#  paste(missing_packages, collapse = ", ")
#  )
#}
#
#message("All required packages available...")
#
#REF_GENOME_FILE <- list.files(
#  path = FILE_GENOME_DIRECTORY,
#  pattern = FILE_GENOME_PATTERN,
#  full.names = TRUE,
#  recursive = TRUE
#)[1]
#
## @QUES: Need to rename this?
#FEATURE_FILE <- list.files(
#  path = FILE_FEATURE_DIRECTORY,
#  pattern = FILE_FEATURE_PATTERN,
#  full.names = TRUE
#)[1]
#
#if (length(FEATURE_FILE) == 0) {
#  warning(sprintf("No feature files found matching pattern '%s' in: %s",
#  FILE_FEATURE_PATTERN,
#  FILE_FEATURE_DIRECTORY
#  ))
#}
#
#if (length(REF_GENOME_FILE) == 0) {
#  stop(sprintf("No reference genome files found matching pattern '%s' in: %s",
#  FILE_GENOME_DIRECTORY,
#  FILE_FEATURE_DIRECTORY
#  ))
#}
#
#REFERENCE_GENOME_DSS <- Biostrings::readDNAStringSet(REF_GENOME_FILE)
#CHROMOSOME_WIDTHS <- REFERENCE_GENOME_DSS[CHROMOSOME_TO_PLOT]@ranges@width
#CHROMOSOMES_IN_ROMAN <- paste0("chr", utils::as.roman(CHROMOSOME_TO_PLOT))
#
## @FIX: Better name.
#ALL_CHROMOSOMES <- REFERENCE_GENOME_DSS@ranges@NAMES
#
#for (chromosome in ALL_CHROMOSOMES) {
#  message("Current chromosome: ", chromosome)
#}
##for (chromosome in ) {
##  message("Current chromosome: ", chromosome)
##
##}
#message("REFERENCE chromosome is same as sample chromosomes: ", as.character(identical(ALL_CHROMOSOMES, names(CHROMOSOME_LENGTHS_nls))))
#
##reads_plus <- READS_CHRI_galn[strand(READS_CHRI_galn) == "+"]
##coverage_plus <- coverage(reads_plus, width = unname(chromosome_length))
##summarizeOverlaps(features = WINDOW_TILES_vct, reads = REFERENCE_BAM_PATH, ignore.strand = TRUE)
