# Attempt to replicate eaton peak calling protocol found in supplementary.
# Eaton et al 2010. CitationKey: Eaton2010conserved
library(MASS)
library(GenomicAlignments)
library(rtracklayer)
library(Rsamtools)

# INPUTS
# REQUIRED
# CHIP reads as bam (turn into grange or galignments)
REFERENCE_BAM_PATH <- "/home/luised94/data/100303Bel/alignment/consolidated_034475_sequence_to_S288C_sorted.bam"
# ORF annotations - filter the gff after loading
# @FIX: Score column has CDS. Need to fix during the import
#orf_bed_file_path <- "~/data/feature_files/20250423_orf_sgd.bed"
GFF_FILE_PATH <- "~/data/feature_files/240830_saccharomyces_cerevisiae.gff"

# OPTIONAL
# INPUT_BAM <- ""
# BLACKLIST_BED_FILE_PATH <- ""

# Configuration
THRESHOLD_FOLD <- 3
CROSS_CORR_MIN_OFFSET_bp <- 100
CROSS_CORR_MAX_OFFSET_bp <- 300
CROSS_CORR_STEP_bp <- 10
CROSS_CORR_BIN_SIZE_bp <- 20
SLIDING_WINDOWS_STEP_SIZE_bp <- 25
PVALUE_THRESHOLD <- 10^-6
WINDOW_MAX_MERGE_GAP_bp <- 50
MIN_REPLICATE_OVERLAP_pct <- 0.5
CANDIDATE_WINDOW_SIZE_bp <- 300
ENRICHMENT_FOLD <- 3
CANDIDATE_REGION_MIN_SIZE_bp <- 500
OFFSETS_bp <- seq(
  from = CROSS_CORR_MIN_OFFSET_bp,
  to = CROSS_CORR_MAX_OFFSET_bp,
  by = CROSS_CORR_STEP_bp
)


# Import the GFF file. rtracklayer will parse all the rich metadata.
# granges
gff_annotations <- import(GFF_FILE_PATH)
orf_annotations <- gff_annotations[mcols(gff_annotations)$type == "CDS"]
REFERENCE_BAM_HEADER <- Rsamtools::scanBamHeader(REFERENCE_BAM_PATH)
CHROMOSOME_LENGTHS_nls <- REFERENCE_BAM_HEADER[[REFERENCE_BAM_PATH]]$targets
CHROMOSOME_NAMES_chr <- names(CHROMOSOME_LENGTHS_nls)
ALL_READS_galn <- GenomicAlignments::readGAlignments(
  REFERENCE_BAM_PATH,
  param = ScanBamParam(
    which = GRanges(
      seqnames = CHROMOSOME_NAMES_chr,
      ranges = IRanges(1, CHROMOSOME_LENGTHS_nls)
    )
  )
)

WINDOW_TILES_vct <- GenomicRanges::tileGenome(
  seqlengths = CHROMOSOME_LENGTHS_nls,
  tilewidth = CANDIDATE_WINDOW_SIZE_bp,
  cut.last.tile.in.chrom = TRUE
)

reads_by_strand <- split(
  ALL_READS_galn,
  strand(ALL_READS_galn)
)
plus_strand_reads <- reads_by_strand$`+`
minus_strand_reads <- reads_by_strand$`-`
plus_gr <- granges(plus_strand_reads)
minus_gr <- granges(minus_strand_reads)
coverage_plus <- coverage(plus_strand_reads, width = unname(CHROMOSOME_LENGTHS_nls))
coverage_minus <- coverage(minus_strand_reads, width = unname(CHROMOSOME_LENGTHS_nls))

se <- GenomicAlignments::summarizeOverlaps(
  features = WINDOW_TILES_vct,
  reads = ALL_READS_galn,
  ignore.strand = TRUE
)
score_by_window <- assay(se)[,1]
mean_score <- mean(score_by_window)
threshold_score <- mean_score * THRESHOLD_FOLD
is_candidate_window <- score_by_window >= threshold_score
scores_above_threshold <- score_by_window[is_candidate_window]
candidate_windows <- WINDOW_TILES_vct[is_candidate_window]
candidate_regions <- reduce(candidate_windows, ignore.strand = TRUE)

overlaps <- GenomicAlignments::findOverlaps(
  query = candidate_regions,
  subject = candidate_windows
)
num_windows <- countQueryHits(overlaps)
mcols(candidate_regions)$num_merged_windows <- num_windows
large_candidate_regions <- candidate_regions[width(candidate_regions) >= CANDIDATE_REGION_MIN_SIZE_bp]

# Create a list to store the results
cross_corr_results <- lapply(seq_along(large_candidate_regions), function(i) {
  region <- large_candidate_regions[i]
  chrom <- as.character(seqnames(region))

  # --- Bin coverage into 20bp bins ---
  bins <- unlist(tile(
    region,
    width = CROSS_CORR_BIN_SIZE_bp
  ))
  ir_bins <- ranges(bins)
  
  # Get the coverage views for the bins on the appropriate chromosome
  plus_views <- Views(coverage_plus[[chrom]], ir_bins)
  minus_views <- Views(coverage_minus[[chrom]], ir_bins)
  
  # Calculate mean coverage in each bin
  plus_binned_cov <- viewMeans(plus_views)
  minus_binned_cov <- viewMeans(minus_views)

  # --- Shifting and Correlation ---
  correlations <- sapply(OFFSETS_bp, function(offset) {
    # Convert bp offset to bin offset. Use floor to ensure integer.
    bin_offset <- floor(offset / 20)
    
    # Shift the minus strand coverage vector upstream (left shift)
    # We need to pad the end with NAs and remove elements from the start
    n_bins <- length(minus_binned_cov)
    if (bin_offset >= n_bins) return(NA) # Offset is larger than the vector
    
    shifted_minus <- c(
      minus_binned_cov[(bin_offset + 1):n_bins],
      rep(NA, bin_offset)
    )
    
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
    best_offset <- OFFSETS_bp[max_idx]
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

# --- Inspection ---
#print("Summary of the imported GFF annotations:")
#print(gff_annotations)

# Inspect the feature types available in the file
#print("Available feature types in the GFF file:")
#print(table(mcols(gff_annotations)$type))

# --- 1. Generate Sliding Windows ---
# Define the step size for the sliding windows
# Create a GRanges object representing all sliding windows across the genome
# This works by creating all start sites and then defining ranges of 'window_size'
window_starts <- lapply(
  names(CHROMOSOME_LENGTHS_nls),
  function(chrom) {
    seq(
      from = 1,
      to = CHROMOSOME_LENGTHS_nls[chrom] - window_size + 1,
      by = SLIDING_WINDOWS_STEP_SIZE_bp
    )
  }
)

# Create proper Seqinfo object with lengths
seqinfo_with_lengths <- Seqinfo(
  seqnames = names(CHROMOSOME_LENGTHS_nls),
  seqlengths = CHROMOSOME_LENGTHS_nls
)

sliding_windows <- GRanges(
  seqnames = rep(
    names(CHROMOSOME_LENGTHS_nls),
    lengths(window_starts)
  ),
  ranges = IRanges(
    start = unlist(window_starts), 
    width = window_size
  ),
  # Use seqinfo from annotations for consistency
  seqinfo = seqinfo_with_lengths
)

# Count centered reads in each sliding window
# fragment_centers are the strand-agnostic shifted read positions
# sliding_windows are genome-wide windows of width w, step 25bp
SLIDING_WINDOW_COUNTS_se <- GenomicAlignments::summarizeOverlaps(
  features = sliding_windows,
  reads = fragment_centers,
  mode = "Union",
  inter.feature = FALSE,
  ignore.strand = TRUE,
  singleEnd = TRUE  # CRITICAL: treat fragment centers as single-end
)

# Extract counts as a simple vector
SLIDING_WINDOW_COUNTS_vct <- assay(SLIDING_WINDOW_COUNTS_se)[, 1]


# Find all overlaps between sliding windows and ORF annotations (CDS)
# This returns a Hits object mapping query (windows) to subject (ORFs)
WINDOW_ORF_OVERLAPS <- GenomicRanges::findOverlaps(
  query = sliding_windows,
  subject = orf_annotations,
  ignore.strand = TRUE
)

# Calculate the actual intersection ranges to measure overlap width
# pintersect gives us the overlapping portions
OVERLAP_RANGES_gr <- GenomicRanges::pintersect(
  sliding_windows[queryHits(WINDOW_ORF_OVERLAPS)],
  orf_annotations[subjectHits(WINDOW_ORF_OVERLAPS)]
)

# Calculate overlap width for each hit
overlap_widths_vct <- width(OVERLAP_RANGES_gr)

# Aggregate total overlap width per window
# Using tapply to sum overlap widths grouped by window index
overlap_sums_by_window <- tapply(
  overlap_widths_vct,
  queryHits(WINDOW_ORF_OVERLAPS),
  sum
)

# Note: tapply returns a named vector with window indices as names
# We need to expand this to all windows (most have 0 overlap)

# Create a vector for all windows, initialize with 0
total_overlap_width_vct <- rep(0, length(sliding_windows))

# Fill in the overlaps we found
window_indices_with_overlap <- as.integer(names(overlap_sums_by_window))
total_overlap_width_vct[window_indices_with_overlap] <- overlap_sums_by_window

# Cap at window_size to handle any potential double-counting
total_overlap_width_vct <- pmin(total_overlap_width_vct, window_size)

# Calculate overlap percentage
overlap_percentage_vct <- total_overlap_width_vct / window_size
# Filter windows: retain only those with ò50% ORF overlap

is_background_window <- overlap_percentage_vct >= 0.5
BACKGROUND_WINDOWS_gr <- sliding_windows[is_background_window]
BACKGROUND_WINDOW_COUNTS_vct <- SLIDING_WINDOW_COUNTS_vct[is_background_window]

valid_counts <- BACKGROUND_WINDOW_COUNTS_vct[!is.na(BACKGROUND_WINDOW_COUNTS_vct)]

nb_fit <- MASS::fitdistr(
  x = valid_counts,
  densfun = "negative binomial"
)

# Extract parameters
NB_MU <- nb_fit$estimate["mu"]
NB_THETA <- nb_fit$estimate["size"] 

mom_mu <- mean(valid_counts)
mom_var <- var(valid_counts)

# For negative binomial: var = æ + æý/?
# Solving for ?: ? = æý / (var - æ)
mom_theta <- mom_mu^2 / (mom_var - mom_mu)


# Load the published peaks (adjust path as needed)
#PAPER_PEAKS_PATH <- "~/data/feature_files/240830_eaton_peaks.bed"
#paper_peaks <- rtracklayer::import(PAPER_PEAKS_PATH, format = "bed")
#
#cat("Our peaks:", length(MERGED_PEAKS_gr), "\n")
#cat("Paper peaks:", length(paper_peaks), "\n")
#
#overlaps <- findOverlaps(MERGED_PEAKS_gr, paper_peaks)
#cat("Our peaks overlapping paper:", length(unique(queryHits(overlaps))), "\n")
#cat("Paper peaks overlapping ours:", length(unique(subjectHits(overlaps))), "\n")
