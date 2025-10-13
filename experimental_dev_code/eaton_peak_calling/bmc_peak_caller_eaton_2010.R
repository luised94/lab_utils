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

# Goodness-of-fit diagnostics
theoretical_mean <- NB_MU
theoretical_var <- NB_MU + (NB_MU^2 / NB_THETA)
# 2. Overdispersion parameter (variance-to-mean ratio)
cat("Overdispersion (var/mean):\n")
cat("  Observed:", var(valid_counts) / mean(valid_counts), "\n")
cat("  Theoretical:", theoretical_var / theoretical_mean, "\n\n")

et.seed(42)
simulated_counts <- rnbinom(
  n = length(valid_counts),
  mu = NB_MU,
  size = NB_THETA
)

cat("Quantile comparison (observed vs theoretical):\n")
quantiles <- c(0.25, 0.50, 0.75, 0.90, 0.95, 0.99)
obs_quantiles <- quantile(valid_counts, probs = quantiles)
theo_quantiles <- quantile(simulated_counts, probs = quantiles)

comparison_df <- data.frame(
  Quantile = quantiles,
  Observed = obs_quantiles,
  Theoretical = theo_quantiles,
  Difference = obs_quantiles - theo_quantiles
)
print(comparison_df)

# STEPS 11-12: GENOME-WIDE SCANNING AND STATISTICAL TESTING
# For each window count k, calculate P(X ò k) where X ~ NegBinom(æ, ?)
# pnbinom with lower.tail=FALSE gives P(X > k), so we use P(X > k-1) = P(X ò k)
WINDOW_PVALUES_vct <- pnbinom(
  q = SLIDING_WINDOW_COUNTS_vct - 1,  # q is the quantile
  mu = NB_MU,                          # mean from fitted model
  size = NB_THETA,                     # dispersion from fitted model
  lower.tail = FALSE                   # upper tail: P(X > q)
)

cat("Summary of p-values:\n")
print(summary(WINDOW_PVALUES_vct))

# Apply significance threshold from configuration
IS_SIGNIFICANT_WINDOW <- WINDOW_PVALUES_vct <= PVALUE_THRESHOLD

# Extract significant windows and their statistics
SIGNIFICANT_WINDOWS_gr <- sliding_windows[IS_SIGNIFICANT_WINDOW]
SIGNIFICANT_WINDOW_COUNTS_vct <- SLIDING_WINDOW_COUNTS_vct[IS_SIGNIFICANT_WINDOW]
SIGNIFICANT_WINDOW_PVALUES_vct <- WINDOW_PVALUES_vct[IS_SIGNIFICANT_WINDOW]

# Add counts and p-values as metadata
mcols(SIGNIFICANT_WINDOWS_gr)$count <- SIGNIFICANT_WINDOW_COUNTS_vct
mcols(SIGNIFICANT_WINDOWS_gr)$pvalue <- SIGNIFICANT_WINDOW_PVALUES_vct
mcols(SIGNIFICANT_WINDOWS_gr)$neg_log10_pvalue <- -log10(SIGNIFICANT_WINDOW_PVALUES_vct)

# Summary statistics
cat("\n=== Significant Window Summary ===\n")
cat("Total windows tested:", length(sliding_windows), "\n")
cat("Significant windows (p ó", PVALUE_THRESHOLD, "):", length(SIGNIFICANT_WINDOWS_gr), "\n")
cat("Percentage significant:", 
    round(100 * length(SIGNIFICANT_WINDOWS_gr) / length(sliding_windows), 2), "%\n")

cat("\nSignificant window counts:\n")
print(summary(SIGNIFICANT_WINDOW_COUNTS_vct))

cat("\nDistribution across chromosomes:\n")
print(table(seqnames(SIGNIFICANT_WINDOWS_gr)))

cat("\nSteps 11-12 complete!\n")

# STEP 13: MERGE SIGNIFICANT WINDOWS INTO PEAKS
# Reduce merges overlapping/nearby ranges
# min.gapwidth controls maximum gap: gaps < min.gapwidth are merged
# Setting min.gapwidth = WINDOW_MAX_MERGE_GAP_bp + 1 means:
# "merge windows if gap between them is ó WINDOW_MAX_MERGE_GAP_bp"
MERGED_PEAKS_gr <- GenomicRanges::reduce(
  SIGNIFICANT_WINDOWS_gr,
  min.gapwidth = WINDOW_MAX_MERGE_GAP_bp + 1,
  ignore.strand = TRUE
)

# For each merged peak, calculate summary statistics
# We need to find which original windows contributed to each peak
peak_overlaps <- GenomicRanges::findOverlaps(
  query = MERGED_PEAKS_gr,
  subject = SIGNIFICANT_WINDOWS_gr
)
# For each peak, aggregate statistics from contributing windows
peak_stats <- lapply(seq_along(MERGED_PEAKS_gr), function(i) {
  # Get indices of windows that overlap this peak
  window_indices <- subjectHits(peak_overlaps)[queryHits(peak_overlaps) == i]
  
  # Get the corresponding counts and p-values
  counts <- SIGNIFICANT_WINDOW_COUNTS_vct[window_indices]
  pvalues <- SIGNIFICANT_WINDOW_PVALUES_vct[window_indices]
  neg_log10_pvals <- -log10(pvalues)
  
  # Return summary statistics
  data.frame(
    max_count = max(counts),
    mean_count = mean(counts),
    max_neg_log10_pvalue = max(neg_log10_pvals),
    num_windows = length(window_indices)
  )
})

# Convert to DataFrame and add as metadata
peak_stats_df <- do.call(rbind, peak_stats)
mcols(MERGED_PEAKS_gr) <- DataFrame(peak_stats_df)

# Add score column (peak score = maximum -log10 p-value)
mcols(MERGED_PEAKS_gr)$score <- mcols(MERGED_PEAKS_gr)$max_neg_log10_pvalue

cat("=== Peak Summary Statistics ===\n")
cat("Peak width distribution:\n")
print(summary(width(MERGED_PEAKS_gr)))

cat("\nPeak score distribution (max -log10 p-value):\n")
print(summary(mcols(MERGED_PEAKS_gr)$score))

cat("\nWindows per peak:\n")
print(summary(mcols(MERGED_PEAKS_gr)$num_windows))

cat("\nPeaks per chromosome:\n")
print(table(seqnames(MERGED_PEAKS_gr)))

cat("\nStep 13 complete!\n")

# Load the published peaks (adjust path as needed)
# Convert paper peak seqnames to match ours (Roman numerals with "chr" prefix)
# Mapping: 1chrI, 2chrII, etc.
roman_numerals <- c("I", "II", "III", "IV", "V", "VI", "VII", "VIII", 
                    "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI")
# Create mapping from numbers to chrRoman
seqname_map <- setNames(
  paste0("chr", roman_numerals),
  as.character(1:16)
)
PAPER_PEAKS_PATH <- "~/data/feature_files/240830_eaton_peaks.bed"
paper_peaks <- rtracklayer::import(PAPER_PEAKS_PATH, format = "bed")

cat("Our peaks:", length(MERGED_PEAKS_gr), "\n")
cat("Paper peaks:", length(paper_peaks), "\n")
cat("Our peaks overlapping paper:", length(unique(queryHits(overlaps))), "\n")
cat("Paper peaks overlapping ours:", length(unique(subjectHits(overlaps))), "\n")

# Apply mapping to paper peaks
paper_peaks_renamed <- paper_peaks
seqlevels(paper_peaks_renamed) <- seqname_map[seqlevels(paper_peaks_renamed)]

# Now compare
overlaps <- findOverlaps(MERGED_PEAKS_gr, paper_peaks_renamed)

cat("=== Peak Comparison: Our Implementation vs Paper ===\n")
cat("Our peaks:", length(MERGED_PEAKS_gr), "\n")
cat("Paper peaks:", length(paper_peaks_renamed), "\n")
cat("Our peaks overlapping paper:", length(unique(queryHits(overlaps))), "\n")
cat("Paper peaks overlapping ours:", length(unique(subjectHits(overlaps))), "\n")

# Calculate overlap percentages
our_overlap_pct <- 100 * length(unique(queryHits(overlaps))) / length(MERGED_PEAKS_gr)
paper_overlap_pct <- 100 * length(unique(subjectHits(overlaps))) / length(paper_peaks_renamed)

cat("\nOverlap rate:\n")
cat("  ", round(our_overlap_pct, 1), "% of our peaks overlap paper\n")
cat("  ", round(paper_overlap_pct, 1), "% of paper peaks overlap ours\n")

# Peaks unique to each set
our_unique <- length(MERGED_PEAKS_gr) - length(unique(queryHits(overlaps)))
paper_unique <- length(paper_peaks_renamed) - length(unique(subjectHits(overlaps)))

cat("\nUnique peaks:\n")
cat("  Our unique peaks:", our_unique, "\n")
cat("  Paper unique peaks:", paper_unique, "\n")
# THIS SECTION IS UNTESTED. USE FOR REFERENCE.
# ============================================================================
# STEPS 14-15: REPLICATE CONSENSUS (for separate replicates)
# ============================================================================

# Assume you've run the pipeline on 3 separate replicates and have:
# replicate1_peaks_gr, replicate2_peaks_gr, replicate3_peaks_gr

# Put all replicate peak sets in a list
replicate_peaks_list <- list(
  rep1 = replicate1_peaks_gr,
  rep2 = replicate2_peaks_gr,
  rep3 = replicate3_peaks_gr
)

# --- STEP 14: Find Overlapping Peaks Between Replicates ---

# Function to calculate reciprocal overlap between two peak sets
calculate_reciprocal_overlap <- function(peaks_query, peaks_subject) {
  overlaps <- findOverlaps(peaks_query, peaks_subject)
  
  # For each overlap, calculate reciprocal overlap percentage
  reciprocal_overlaps <- sapply(seq_along(overlaps), function(i) {
    query_idx <- queryHits(overlaps)[i]
    subject_idx <- subjectHits(overlaps)[i]
    
    # Get the two peaks
    query_peak <- peaks_query[query_idx]
    subject_peak <- peaks_subject[subject_idx]
    
    # Calculate intersection
    intersection <- intersect(query_peak, subject_peak)
    overlap_width <- sum(width(intersection))
    
    # Reciprocal overlap = overlap / min(length1, length2)
    min_width <- min(width(query_peak), width(subject_peak))
    reciprocal_overlap <- overlap_width / min_width
    
    return(reciprocal_overlap)
  })
  
  # Keep only overlaps meeting threshold
  passing_overlaps <- overlaps[reciprocal_overlaps >= MIN_REPLICATE_OVERLAP_pct]
  
  return(passing_overlaps)
}
# --- STEP 15: Build Consensus Peak Set ---

# Strategy depends on number of replicates:
n_replicates <- length(replicate_peaks_list)

if (n_replicates == 2) {
  # With 2 replicates: require both
  cat("Finding consensus peaks between 2 replicates...\n")
  
  overlaps_1_2 <- calculate_reciprocal_overlap(
    replicate_peaks_list[[1]], 
    replicate_peaks_list[[2]]
  )
  
  # Peaks from rep1 that overlap rep2
  consensus_indices_rep1 <- unique(queryHits(overlaps_1_2))
  consensus_peaks_rep1 <- replicate_peaks_list[[1]][consensus_indices_rep1]
  
  # For each consensus peak, take the union of coordinates from both replicates
  consensus_peaks_list <- lapply(consensus_indices_rep1, function(idx) {
    # Find all rep2 peaks that overlap this rep1 peak
    overlapping_rep2_indices <- subjectHits(overlaps_1_2)[queryHits(overlaps_1_2) == idx]
    
    # Union of coordinates
    all_peaks <- c(
      replicate_peaks_list[[1]][idx],
      replicate_peaks_list[[2]][overlapping_rep2_indices]
    )
    reduced_peak <- reduce(all_peaks)
    
    # Average the scores
    avg_score <- mean(c(
      mcols(replicate_peaks_list[[1]][idx])$score,
      mcols(replicate_peaks_list[[2]][overlapping_rep2_indices])$score
    ))
    
    mcols(reduced_peak)$score <- avg_score
    mcols(reduced_peak)$num_replicates <- 2
    
    return(reduced_peak)
  })
  
  CONSENSUS_PEAKS_gr <- do.call(c, consensus_peaks_list)
  
} else if (n_replicates >= 3) {
  # With 3+ replicates: require presence in at least 2 (majority)
  cat("Finding consensus peaks across", n_replicates, "replicates...\n")
  cat("Requiring presence in >= 2 replicates\n")
  
  # Create all pairwise comparisons
  # For each peak in each replicate, count how many other replicates it overlaps
  
  # Flatten all peaks with replicate ID
  all_peaks_with_id <- lapply(seq_along(replicate_peaks_list), function(i) {
    peaks <- replicate_peaks_list[[i]]
    mcols(peaks)$replicate_id <- i
    mcols(peaks)$original_index <- seq_along(peaks)
    return(peaks)
  })
  all_peaks_combined <- do.call(c, all_peaks_with_id)
  
  # Find all reciprocal overlaps
  self_overlaps <- findOverlaps(all_peaks_combined, all_peaks_combined)
  
  # Remove self-hits and filter by reciprocal overlap threshold
  # (implementation details for reciprocal overlap calculation...)
  # This is complex - simplified version below
  
  # Simple approach: use reduce to merge all peaks, then count support
  all_peaks_merged <- reduce(all_peaks_combined)
  
  # For each merged peak, count how many replicates contributed
  peak_support <- sapply(seq_along(all_peaks_merged), function(i) {
    merged_peak <- all_peaks_merged[i]
    
    # Count unique replicates overlapping this region
    replicates_overlapping <- sapply(replicate_peaks_list, function(rep_peaks) {
      any(overlapsAny(merged_peak, rep_peaks, minoverlap = 
        as.integer(MIN_REPLICATE_OVERLAP_pct * width(merged_peak))))
    })
    
    sum(replicates_overlapping)
  })
  
  # Keep peaks with support from >= 2 replicates
  CONSENSUS_PEAKS_gr <- all_peaks_merged[peak_support >= 2]
  mcols(CONSENSUS_PEAKS_gr)$num_replicates <- peak_support[peak_support >= 2]
}

cat("\n=== Consensus Peak Summary ===\n")
cat("Total consensus peaks:", length(CONSENSUS_PEAKS_gr), "\n")
cat("Original peaks per replicate:\n")
for (i in seq_along(replicate_peaks_list)) {
  cat("  Replicate", i, ":", length(replicate_peaks_list[[i]]), "peaks\n")
}

