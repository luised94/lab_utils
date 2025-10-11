library(GenomicAlignments)
library(rtracklayer)
library(Rsamtools)

REFERENCE_BAM <- "/home/luised94/data/100303Bel/alignment/consolidated_034475_sequence_to_S288C_sorted.bam"
reference_bam_header <- Rsamtools::scanBamHeader(REFERENCE_BAM)
sample_chromosomes <- names(reference_bam_header[[REFERENCE_BAM]]$targets)
chromosome_name <- sample_chromosomes[1]
chromosome_length <- chromosome_lengths[1]
which_chr <- GRanges(
  chromosome_name, 
  IRanges(1, reference_bam_header[[REFERENCE_BAM]]$targets[[chromosome_name]])
)

reads_chrI <- readGAlignments(
    REFERENCE_BAM, 
    param = ScanBamParam(which = which_chr)
)
reads_by_strand <- split(reads_chrI, strand(reads_chrI))
coverage_plus <- coverage(reads_by_strand$`+`, width = unname(chromosome_length))
coverage_minus <- coverage(reads_by_strand$`-`, width = unname(chromosome_length))
windows_300bp <- tileGenome(
  chromosome_length, 
  tilewidth = 300,
  cut.last.tile.in.chrom = TRUE
)
se <- summarizeOverlaps(
  features = windows_300bp, 
  reads = reads_chrI,
  ignore.strand = TRUE
)
mean_score <- mean(assay(se)[,1])
threshold_fold <- 3
threshold_score <- mean_score * threshold_fold
is_candidate_window <- assay(se)[,1] >= threshold_score
scores_above_threshold <- assay(se)[,1][is_candidate_window]
candidate_windows <- windows_300bp[is_candidate_window]
candidate_regions <- reduce(candidate_windows, ignore.strand = TRUE)

overlaps <- findOverlaps(candidate_regions, candidate_windows)
num_windows <- countQueryHits(overlaps)
mcols(candidate_regions)$num_merged_windows <- num_windows
large_candidate_regions <- candidate_regions[width(candidate_regions) >= 500]

# Define offsets in base pairs
offsets <- seq(100, 300, by = 10)

# Create a list to store the results
cross_corr_results <- lapply(seq_along(large_candidate_regions), function(i) {
    region <- large_candidate_regions[i]
    chrom <- as.character(seqnames(region))

    # --- Bin coverage into 20bp bins ---
    bins <- unlist(tile(region, width = 20))
    ir_bins <- ranges(bins)
    
    # Get the coverage views for the bins on the appropriate chromosome
    plus_views <- Views(coverage_plus[[chrom]], ir_bins)
    minus_views <- Views(coverage_minus[[chrom]], ir_bins)
    
    # Calculate mean coverage in each bin
    plus_binned_cov <- viewMeans(plus_views)
    minus_binned_cov <- viewMeans(minus_views)

    # --- Shifting and Correlation ---
    correlations <- sapply(offsets, function(offset) {
        # Convert bp offset to bin offset. Use floor to ensure integer.
        bin_offset <- floor(offset / 20)
        
        # Shift the minus strand coverage vector upstream (left shift)
        # We need to pad the end with NAs and remove elements from the start
        n_bins <- length(minus_binned_cov)
        if (bin_offset >= n_bins) return(NA) # Offset is larger than the vector
        
        shifted_minus <- c(minus_binned_cov[(bin_offset + 1):n_bins], rep(NA, bin_offset))
        
        # Calculate Pearson correlation, removing NAs
        cor_test <- cor.test(plus_binned_cov, shifted_minus, method = "pearson")
        return(cor_test$estimate)
    })

    # Find the offset with the maximum correlation
    if (all(is.na(correlations))) {
        max_corr <- NA
        best_offset <- NA
    } else {
        max_idx <- which.max(correlations)
        max_corr <- correlations[max_idx]
        best_offset <- offsets[max_idx]
    }
    
    return(data.frame(
        region_id = i, 
        max_correlation = max_corr, 
        best_offset = best_offset
    ))
})

# Combine results into a single data frame
final_results <- do.call(rbind, cross_corr_results)

# You can add these results back to your GRanges object
mcols(large_candidate_regions) <- cbind(mcols(large_candidate_regions), final_results)

median_offset <- median(final_results$best_offset)
fragment_size <- median_offset
window_size <- fragment_size * 2

# 1. Calculate the shift amount
shift_amount <- floor(fragment_size / 2)

# 2. Split reads by strand
reads_by_strand <- split(reads_chrI, strand(reads_chrI))
plus_strand_reads <- reads_by_strand$`+`
minus_strand_reads <- reads_by_strand$`-`
plus_gr <- granges(plus_strand_reads)
minus_gr <- granges(minus_strand_reads)

# 4. Resize and shift, creating potentially out-of-bounds ranges
shifted_plus_raw <- shift(resize(plus_gr, width = 1, fix = "start"), shift = shift_amount)
shifted_minus_raw <- shift(resize(minus_gr, width = 1, fix = "end"), shift = -shift_amount)

trimmed_plus <- trim(shifted_plus_raw)
trimmed_minus <- trim(shifted_minus_raw)
out_of_bounds_plus <- length(shifted_plus_raw) - length(trimmed_plus)
out_of_bounds_minus <- length(shifted_minus_raw) - length(trimmed_minus)
cat("Number of out-of-bounds plus-strand fragments:", out_of_bounds_plus, "\n")
cat("Number of out-of-bounds minus-strand fragments:", out_of_bounds_minus, "\n")
fragment_centers <- c(trimmed_plus, trimmed_minus)
strand(fragment_centers) <- "*"
# Import the BED file into a GRanges object
# Score column has CDS. Need to fix during the import
#orf_bed_file_path <- "~/data/feature_files/20250423_orf_sgd.bed"
#orf_annotations <- import(orf_bed_file_path)
gff_file_path <- "~/data/feature_files/240830_saccharomyces_cerevisiae.gff"

# Import the GFF file. rtracklayer will parse all the rich metadata.
gff_annotations <- import(gff_file_path)

# --- Inspection ---
print("Summary of the imported GFF annotations:")
print(gff_annotations)

# Inspect the feature types available in the file
print("Available feature types in the GFF file:")
print(table(mcols(gff_annotations)$type))

# --- 1. Generate Sliding Windows ---
# Define the step size for the sliding windows
step_size <- 25

# Create a GRanges object representing all sliding windows across the genome
# This works by creating all start sites and then defining ranges of 'window_size'
window_starts <- lapply(
  names(chromosome_lengths), function(chrom) {
  seq(
    from = 1,
    to = chromosome_lengths[chrom] - window_size + 1,
    by = step_size
  )
})

sliding_windows <- GRanges(
  seqnames = rep(
        names(chromosome_lengths),
        lengths(window_starts)
  ),
  ranges = IRanges(start = unlist(window_starts), width = window_size),
  seqinfo = seqinfo(orf_annotations) # Use seqinfo from annotations for consistency
)

######### JUNK ###########
#CHROMOSOME_TO_PLOT <- 10
#FILE_FEATURE_DIRECTORY <- file.path(Sys.getenv("HOME"), "data", "feature_files")
#FILE_GENOME_DIRECTORY <- file.path(Sys.getenv("HOME"), "data", "REFGENS")
#
#FILE_GENOME_PATTERN <- "S288C_refgenome.fna"
#FILE_FEATURE_PATTERN <- "eaton_peaks"
## @ques: maybe I should download this file locally.
#REFERENCE_BAM <- "/home/luised94/data/100303Bel/alignment/consolidated_034475_sequence_to_S288C_sorted.bam"
#reference_bam_header <- Rsamtools::scanBamHeader(REFERENCE_BAM)
#print(reference_bam_header[[REFERENCE_BAM]]$targets)
#
#stopifnot(
#  #"REFERENCE_BAM does not exist." =
#  #file.exists(REFERENCE_BAM),
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
#    "Missing packages. Please install using renv:\n",
#    paste(missing_packages, collapse = ", ")
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
#    FILE_FEATURE_PATTERN,
#    FILE_FEATURE_DIRECTORY
#  ))
#}
#
#if (length(REF_GENOME_FILE) == 0) {
#  stop(sprintf("No reference genome files found matching pattern '%s' in: %s",
#    FILE_GENOME_DIRECTORY,
#    FILE_FEATURE_DIRECTORY
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
#message("REFERENCE chromosome is same as sample chromosomes: ", as.character(identical(ALL_CHROMOSOMES, names(reference_bam_header[[REFERENCE_BAM]]$targets))))
#chromosome_lengths <- reference_bam_header[[REFERENCE_BAM]]$targets
#
##reads_plus <- reads_chrI[strand(reads_chrI) == "+"]
##coverage_plus <- coverage(reads_plus, width = unname(chromosome_length))
##summarizeOverlaps(features = windows_300bp, reads = REFERENCE_BAM, ignore.strand = TRUE)
