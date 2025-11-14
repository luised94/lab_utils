# ============================================================================
# SCRIPT: CRISPR Guide RNA Quality Control Analysis
# PURPOSE: Add quality metrics to guide RNA dataset
# INPUT: grna_results.tsv (from sgd_ngg_analysis.R pipeline)
# OUTPUT: Full dataset + top 5 per gene + summary statistics
# ============================================================================

# ============================================================================
# SECTION: Configuration constants
# ============================================================================

HOME_DIRECTORY <- Sys.getenv("HOME")
OUTPUT_DIRECTORY_NAME <- "data/fasta_files"
REFERENCE_GENOME_DIR_NAME <- "data/reference_genomes"

INPUT_FILE_NAME_chr <- "grna_results.tsv"
GENOME_FILE_NAME_chr <- "Saccharomyces_cerevisiae_W303_GCA0021635151_genome.fna"
OUTPUT_FULL_NAME_chr <- "grna_results_quality.tsv"
OUTPUT_TOP5_NAME_chr <- "grna_top5_per_gene.tsv"
OUTPUT_SUMMARY_NAME_chr <- "grna_quality_summary.txt"

# Directory paths
OUTPUT_DIRECTORY_PATH <- file.path(HOME_DIRECTORY, OUTPUT_DIRECTORY_NAME)
INPUT_DIRECTORY_PATH <- file.path(HOME_DIRECTORY, OUTPUT_DIRECTORY_NAME)
REFERENCE_GENOME_DIR_PATH <- file.path(HOME_DIRECTORY, REFERENCE_GENOME_DIR_NAME)

# File paths
INPUT_FILE_PATH_chr <- file.path(OUTPUT_DIRECTORY_PATH, "grna_results.tsv")
GENOME_FILE_PATH_chr <- file.path(REFERENCE_GENOME_DIR_PATH, "Saccharomyces_cerevisiae_W303_GCA0021635151_genome.fna")
OUTPUT_FULL_PATH_chr <- file.path(OUTPUT_DIRECTORY_PATH, "grna_results_quality.tsv")
OUTPUT_TOP5_PATH_chr <- file.path(OUTPUT_DIRECTORY_PATH, "grna_top5_per_gene.tsv")
OUTPUT_SUMMARY_PATH_chr <- file.path(OUTPUT_DIRECTORY_PATH,"grna_quality_summary.txt")

# GC content thresholds (percent)
GC_MIN_PERCENT_dbl <- 40.0
GC_MAX_PERCENT_dbl <- 60.0
GC_OPTIMAL_PERCENT_dbl <- 50.0  # Target for ranking

# Homopolymer thresholds
HOMOPOLYMER_THRESHOLD_bp_int <- 4L  # Flag if >= 4 consecutive nucleotides
POLYT_THRESHOLD_bp_int <- 4L        # Flag if >= 4 consecutive T's

# Ranking parameters
TOP_N_GUIDES_int <- 5L  # Top guides to output per gene

# Verification
cat("=== Configuration Loaded ===\n")
cat("Input file:", INPUT_FILE_PATH_chr, "\n")
cat("Genome file:", GENOME_FILE_PATH_chr, "\n")
cat("GC range:", GC_MIN_PERCENT_dbl, "-", GC_MAX_PERCENT_dbl, "%\n")
cat("Homopolymer threshold:", HOMOPOLYMER_THRESHOLD_bp_int, "bp\n")
cat("Top guides per gene:", TOP_N_GUIDES_int, "\n")
cat("\n")
# ============================================================================
# SECTION: Load input data
# ============================================================================

# Load guide RNA data from pipeline output
guides_df <- read.table(
  file = INPUT_FILE_PATH_chr,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  na.strings = c("NA", ""),
  comment.char = "",
  quote = ""
)

# Verify data structure
cat("=== Data Loading Summary ===\n")
cat("Total guides loaded:", nrow(guides_df), "\n")
cat("Total columns:", ncol(guides_df), "\n")
cat("Unique genes:", length(unique(guides_df$gene_name)), "\n")
cat("\n")

cat("Column names:\n")
print(names(guides_df))
cat("\n")

cat("Data structure:\n")
str(guides_df, give.attr = FALSE)
cat("\n")

cat("First few rows:\n")
print(head(guides_df, n = 3))
cat("\n")

cat("Strand distribution:\n")
print(table(guides_df$guide_strand))
cat("\n")

cat("CDS overlap distribution:\n")
print(table(guides_df$overlaps_cds, useNA = "ifany"))
cat("\n")

cat("Silencing difficulty distribution:\n")
print(table(guides_df$silencing_difficulty, useNA = "ifany"))
cat("\n")

# ============================================================================
# SECTION: Calculate GC content
# ============================================================================

# Extract guide sequences
guide_sequences_chr <- guides_df$guide_sequence

# Count G and C nucleotides in each guide
gc_count_int <- vapply(
  X = guide_sequences_chr,
  FUN = function(seq_chr) {
    # Count G's and C's
    g_count <- lengths(regmatches(x = seq_chr, m = gregexpr(pattern = "G", text = seq_chr)))
    c_count <- lengths(regmatches(x = seq_chr, m = gregexpr(pattern = "C", text = seq_chr)))
    return(g_count + c_count)
  },
  FUN.VALUE = integer(1),
  USE.NAMES = FALSE
)

# Calculate total length of each guide (should all be 20bp, but verify)
guide_length_int <- nchar(guide_sequences_chr)

# Calculate GC percentage
gc_percent_dbl <- (gc_count_int / guide_length_int) * 100

# Calculate distance from optimal GC content
gc_distance_from_optimal_dbl <- abs(gc_percent_dbl - GC_OPTIMAL_PERCENT_dbl)

# Flag guides within optimal range
gc_in_optimal_range_lgl <- (gc_percent_dbl >= GC_MIN_PERCENT_dbl) & 
                            (gc_percent_dbl <= GC_MAX_PERCENT_dbl)

# Add to dataframe
guides_df$gc_percent_dbl <- gc_percent_dbl
guides_df$gc_distance_from_optimal_dbl <- gc_distance_from_optimal_dbl
guides_df$gc_in_optimal_range_lgl <- gc_in_optimal_range_lgl

# Verify calculations
cat("=== GC Content Analysis ===\n")
cat("Guide length range:", range(guide_length_int), "bp\n")
cat("GC% range:", round(range(gc_percent_dbl), 1), "\n")
cat("GC% mean:", round(mean(gc_percent_dbl), 1), "\n")
cat("GC% median:", round(median(gc_percent_dbl), 1), "\n")
cat("\n")

cat("GC% distribution:\n")
print(summary(gc_percent_dbl))
cat("\n")

cat("Guides in optimal range (40-60%):", sum(gc_in_optimal_range_lgl), 
    "/", nrow(guides_df), 
    "(", round(100 * mean(gc_in_optimal_range_lgl), 1), "%)\n")
cat("Guides below 40%:", sum(gc_percent_dbl < GC_MIN_PERCENT_dbl), "\n")
cat("Guides above 60%:", sum(gc_percent_dbl > GC_MAX_PERCENT_dbl), "\n")
cat("\n")

cat("Sample guides with GC content:\n")
print(guides_df[1:5, c("gene_name", "guide_sequence", "gc_percent_dbl", "gc_in_optimal_range_lgl")])
cat("\n")
# ============================================================================
# SECTION: Detect poly-T runs
# ============================================================================

# Function to find longest homopolymer run of a specific nucleotide
find_longest_run <- function(sequence_chr, nucleotide_chr) {
  # Create pattern for consecutive runs of the nucleotide
  pattern_chr <- paste0(nucleotide_chr, "+")
  
  # Find all runs
  runs_lst <- gregexpr(pattern = pattern_chr, text = sequence_chr)
  run_lengths_int <- attr(runs_lst[[1]], "match.length")
  
  # Return longest run length (0 if none found)
  if (run_lengths_int[1] == -1) {
    return(0L)
  } else {
    return(max(run_lengths_int))
  }
}

# Calculate longest poly-T run for each guide
longest_polyT_run_int <- vapply(
  X = guide_sequences_chr,
  FUN = function(seq_chr) find_longest_run(sequence_chr = seq_chr, nucleotide_chr = "T"),
  FUN.VALUE = integer(1),
  USE.NAMES = FALSE
)

# Flag guides with problematic poly-T runs
has_polyT_problem_lgl <- longest_polyT_run_int >= POLYT_THRESHOLD_bp_int

# Inverse flag for quality scoring (TRUE = good, no poly-T problem)
no_polyT_lgl <- !has_polyT_problem_lgl

# Add to dataframe
guides_df$longest_polyT_run_int <- longest_polyT_run_int
guides_df$has_polyT_problem_lgl <- has_polyT_problem_lgl
guides_df$no_polyT_lgl <- no_polyT_lgl

# Verify calculations
cat("=== Poly-T Analysis ===\n")
cat("Longest poly-T run range:", range(longest_polyT_run_int), "bp\n")
cat("Guides with poly-T >= 4:", sum(has_polyT_problem_lgl), 
    "/", nrow(guides_df),
    "(", round(100 * mean(has_polyT_problem_lgl), 1), "%)\n")
cat("\n")

cat("Poly-T run distribution:\n")
print(table(longest_polyT_run_int))
cat("\n")

cat("Sample guides with poly-T runs >= 4:\n")
polyT_examples_df <- guides_df[has_polyT_problem_lgl, 
                               c("gene_name", "guide_sequence", "longest_polyT_run_int")]
if (nrow(polyT_examples_df) > 0) {
  print(head(polyT_examples_df, n = 10))
} else {
  cat("No guides with poly-T runs >= 4\n")
}
cat("\n")
# ============================================================================
# SECTION: Detect all homopolymer runs
# ============================================================================

# Calculate longest homopolymer run (any nucleotide) for each guide
longest_homopolymer_run_int <- vapply(
  X = guide_sequences_chr,
  FUN = function(seq_chr) {
    # Check all four nucleotides
    runs_int <- c(
      find_longest_run(sequence_chr = seq_chr, nucleotide_chr = "A"),
      find_longest_run(sequence_chr = seq_chr, nucleotide_chr = "T"),
      find_longest_run(sequence_chr = seq_chr, nucleotide_chr = "G"),
      find_longest_run(sequence_chr = seq_chr, nucleotide_chr = "C")
    )
    return(max(runs_int))
  },
  FUN.VALUE = integer(1),
  USE.NAMES = FALSE
)

# Identify which nucleotide forms the longest run
homopolymer_nucleotide_chr <- vapply(
  X = guide_sequences_chr,
  FUN = function(seq_chr) {
    runs_int <- c(
      A = find_longest_run(sequence_chr = seq_chr, nucleotide_chr = "A"),
      T = find_longest_run(sequence_chr = seq_chr, nucleotide_chr = "T"),
      G = find_longest_run(sequence_chr = seq_chr, nucleotide_chr = "G"),
      C = find_longest_run(sequence_chr = seq_chr, nucleotide_chr = "C")
    )
    max_run_int <- max(runs_int)
    if (max_run_int == 0) {
      return(NA_character_)
    } else {
      return(names(runs_int)[which.max(runs_int)])
    }
  },
  FUN.VALUE = character(1),
  USE.NAMES = FALSE
)

# Flag guides with problematic homopolymer runs
has_homopolymer_problem_lgl <- longest_homopolymer_run_int >= HOMOPOLYMER_THRESHOLD_bp_int

# Inverse flag for quality scoring (TRUE = good, no homopolymer problem)
no_homopolymer_lgl <- !has_homopolymer_problem_lgl

# Add to dataframe
guides_df$longest_homopolymer_run_int <- longest_homopolymer_run_int
guides_df$homopolymer_nucleotide_chr <- homopolymer_nucleotide_chr
guides_df$has_homopolymer_problem_lgl <- has_homopolymer_problem_lgl
guides_df$no_homopolymer_lgl <- no_homopolymer_lgl

# Verify calculations
cat("=== Homopolymer Analysis ===\n")
cat("Longest homopolymer run range:", range(longest_homopolymer_run_int), "bp\n")
cat("Guides with homopolymer >= 4:", sum(has_homopolymer_problem_lgl),
    "/", nrow(guides_df),
    "(", round(100 * mean(has_homopolymer_problem_lgl), 1), "%)\n")
cat("\n")

cat("Homopolymer run distribution:\n")
print(table(longest_homopolymer_run_int))
cat("\n")

cat("Homopolymer nucleotide distribution:\n")
print(table(homopolymer_nucleotide_chr, useNA = "ifany"))
cat("\n")

cat("Sample guides with homopolymer runs >= 4:\n")
homopoly_examples_df <- guides_df[has_homopolymer_problem_lgl,
                                  c("gene_name", "guide_sequence", 
                                    "longest_homopolymer_run_int", 
                                    "homopolymer_nucleotide_chr")]
if (nrow(homopoly_examples_df) > 0) {
  print(head(homopoly_examples_df, n = 10))
} else {
  cat("No guides with homopolymer runs >= 4\n")
}
cat("\n")
# ============================================================================
# SECTION: Composite quality scoring
# ============================================================================

# Count quality flags passed (0-3)
quality_flags_passed_int <- (
  as.integer(guides_df$gc_in_optimal_range_lgl) +
  as.integer(guides_df$no_polyT_lgl) +
  as.integer(guides_df$no_homopolymer_lgl)
)

# Add to dataframe
guides_df$quality_flags_passed_int <- quality_flags_passed_int

# Create numeric silencing difficulty score for sorting
# (lower is better: easy=1, moderate=2, difficult=3, very_difficult=4)
silencing_difficulty_score_int <- vapply(
  X = guides_df$silencing_difficulty,
  FUN = function(diff_chr) {
    switch(
      EXPR = diff_chr,
      "easy" = 1L,
      "moderate" = 2L,
      "difficult" = 3L,
      "very_difficult" = 4L,
      NA_integer_  # For any unexpected values
    )
  },
  FUN.VALUE = integer(1),
  USE.NAMES = FALSE
)

guides_df$silencing_difficulty_score_int <- silencing_difficulty_score_int

# Calculate overall quality score (higher is better)
# Formula: (quality_flags * 10) - (silencing_difficulty_score * 3) - gc_distance_from_optimal
# This weights: quality flags most important, then silencing difficulty, then GC optimality
overall_quality_score_dbl <- (
  (quality_flags_passed_int * 10) -
  (silencing_difficulty_score_int * 3) -
  guides_df$gc_distance_from_optimal_dbl
)

guides_df$overall_quality_score_dbl <- overall_quality_score_dbl

# Verify scoring
cat("=== Quality Scoring Summary ===\n")
cat("Quality flags passed distribution:\n")
print(table(quality_flags_passed_int))
cat("\n")

cat("Silencing difficulty score distribution:\n")
print(table(silencing_difficulty_score_int, guides_df$silencing_difficulty, useNA = "ifany"))
cat("\n")

cat("Overall quality score range:", round(range(overall_quality_score_dbl), 2), "\n")
cat("Overall quality score mean:", round(mean(overall_quality_score_dbl), 2), "\n")
cat("Overall quality score median:", round(median(overall_quality_score_dbl), 2), "\n")
cat("\n")

cat("Distribution of guides by combined quality:\n")
cat("Perfect quality (all 3 flags + easy silencing):", 
    sum(quality_flags_passed_int == 3 & silencing_difficulty_score_int == 1), "\n")
cat("High quality (2-3 flags + easy/moderate silencing):",
    sum(quality_flags_passed_int >= 2 & silencing_difficulty_score_int <= 2), "\n")
cat("Moderate quality (1-2 flags):",
    sum(quality_flags_passed_int >= 1 & quality_flags_passed_int < 3), "\n")
cat("Low quality (0 flags or difficult silencing):",
    sum(quality_flags_passed_int == 0 | silencing_difficulty_score_int >= 3), "\n")
cat("\n")

cat("Sample guides with quality scores:\n")
sample_indices_int <- c(1, 100, 500, 1000, 1500)
print(guides_df[sample_indices_int, 
                c("gene_name", "guide_sequence", "gc_percent_dbl",
                  "quality_flags_passed_int", "silencing_difficulty",
                  "overall_quality_score_dbl")])
cat("\n")
# ============================================================================
# SECTION: Sort and rank guides
# ============================================================================

# Sort entire dataset by quality score (descending)
guides_sorted_df <- guides_df[order(guides_df$overall_quality_score_dbl, 
                                    decreasing = TRUE), ]

# Reset row names after sorting
rownames(guides_sorted_df) <- NULL

# Add overall rank
guides_sorted_df$overall_rank_int <- seq_len(nrow(guides_sorted_df))

# Calculate per-gene rank
# Split by gene, add rank within each gene, then recombine
guides_with_gene_rank_df <- do.call(
  what = rbind,
  args = lapply(
    X = split(x = guides_sorted_df, f = guides_sorted_df$gene_name),
    FUN = function(gene_df) {
      # Sort by quality score within gene (already sorted, but explicit)
      gene_df <- gene_df[order(gene_df$overall_quality_score_dbl, 
                               decreasing = TRUE), ]
      
      # Add gene-specific rank
      gene_df$gene_rank_int <- seq_len(nrow(gene_df))
      
      return(gene_df)
    }
  )
)

# Reset row names
rownames(guides_with_gene_rank_df) <- NULL

# Update main dataframe
guides_df <- guides_with_gene_rank_df

# Verify ranking
cat("=== Ranking Summary ===\n")
cat("Total guides:", nrow(guides_df), "\n")
cat("\n")

cat("Top 10 guides overall:\n")
print(guides_df[1:10, c("overall_rank_int", "gene_name", "guide_sequence",
                        "overall_quality_score_dbl", "quality_flags_passed_int",
                        "silencing_difficulty")])
cat("\n")

cat("Guides per gene:\n")
guides_per_gene_int <- table(guides_df$gene_name)
print(guides_per_gene_int)
cat("\n")

cat("Sample of gene-specific rankings (RPO21):\n")
rpo21_guides_df <- guides_df[guides_df$gene_name == "RPO21", ]
print(head(rpo21_guides_df[, c("gene_rank_int", "guide_sequence", 
                                "overall_quality_score_dbl",
                                "quality_flags_passed_int",
                                "silencing_difficulty")], n = 10))
cat("\n")
# ============================================================================
# SECTION: Write output files
# ============================================================================

# ============================================================================
# Output 1: Full dataset with quality metrics
# ============================================================================

write.table(
  x = guides_df,
  file = OUTPUT_FULL_PATH_chr,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  na = "NA"
)

cat("=== Output 1: Full Dataset ===\n")
cat("Written to:", OUTPUT_FULL_PATH_chr, "\n")
cat("Total guides:", nrow(guides_df), "\n")
cat("Total columns:", ncol(guides_df), "\n")
cat("\n")

# ============================================================================
# Output 2: Top N guides per gene
# ============================================================================

# Filter to top N guides per gene
top_guides_df <- guides_df[guides_df$gene_rank_int <= TOP_N_GUIDES_int, ]

# Sort by gene name, then by rank within gene
top_guides_df <- top_guides_df[order(top_guides_df$gene_name, 
                                      top_guides_df$gene_rank_int), ]

# Reset row names
rownames(top_guides_df) <- NULL

write.table(
  x = top_guides_df,
  file = OUTPUT_TOP5_PATH_chr,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  na = "NA"
)

cat("=== Output 2: Top Guides Per Gene ===\n")
cat("Written to:", OUTPUT_TOP5_PATH_chr, "\n")
cat("Total guides:", nrow(top_guides_df), "\n")
cat("Guides per gene:", TOP_N_GUIDES_int, "\n")
cat("\n")

cat("Top guide per gene summary:\n")
top1_per_gene_df <- guides_df[guides_df$gene_rank_int == 1, 
                               c("gene_name", "guide_sequence", 
                                 "overall_quality_score_dbl",
                                 "quality_flags_passed_int",
                                 "gc_percent_dbl",
                                 "silencing_difficulty")]
top1_per_gene_df <- top1_per_gene_df[order(top1_per_gene_df$gene_name), ]
rownames(top1_per_gene_df) <- NULL
print(top1_per_gene_df)
cat("\n")
# ============================================================================
# SECTION: Generate summary statistics report
# ============================================================================

# Open connection to summary file
summary_file_con <- file(OUTPUT_SUMMARY_PATH_chr, open = "w")

# Helper function to write section headers
write_section <- function(title_chr) {
  writeLines(text = paste0("\n", paste(rep("=", 80), collapse = ""), "\n"),
             con = summary_file_con)
  writeLines(text = paste0(title_chr, "\n"), con = summary_file_con)
  writeLines(text = paste0(paste(rep("=", 80), collapse = ""), "\n"),
             con = summary_file_con)
}

# ============================================================================
# Report header
# ============================================================================
writeLines(text = "CRISPR GUIDE RNA QUALITY CONTROL REPORT", 
           con = summary_file_con)
writeLines(text = paste("Generated:", Sys.time()), 
           con = summary_file_con)
writeLines(text = paste("Input file:", INPUT_FILE_PATH_chr), 
           con = summary_file_con)

# ============================================================================
# Overall dataset statistics
# ============================================================================
write_section("DATASET OVERVIEW")

writeLines(text = paste("Total guides analyzed:", nrow(guides_df)), 
           con = summary_file_con)
writeLines(text = paste("Total genes:", length(unique(guides_df$gene_name))), 
           con = summary_file_con)
writeLines(text = paste("Average guides per gene:", 
                       round(nrow(guides_df) / length(unique(guides_df$gene_name)), 1)),
           con = summary_file_con)

writeLines(text = "\nGuides per gene:", con = summary_file_con)
guides_per_gene_tbl <- table(guides_df$gene_name)
for (gene_chr in names(guides_per_gene_tbl)) {
  writeLines(text = paste("  ", gene_chr, ":", guides_per_gene_tbl[gene_chr]),
             con = summary_file_con)
}

# ============================================================================
# Quality metrics summary
# ============================================================================
write_section("QUALITY METRICS SUMMARY")

writeLines(text = paste("Quality flags passed distribution:"), 
           con = summary_file_con)
quality_flags_tbl <- table(guides_df$quality_flags_passed_int)
for (flag_count in names(quality_flags_tbl)) {
  writeLines(text = paste("  ", flag_count, "flags:",
                         quality_flags_tbl[flag_count], "guides",
                         paste0("(", round(100 * quality_flags_tbl[flag_count] / nrow(guides_df), 1), "%)")),
             con = summary_file_con)
}

writeLines(text = "\nSilencing difficulty distribution:", 
           con = summary_file_con)
silencing_tbl <- table(guides_df$silencing_difficulty)
for (diff_chr in c("easy", "moderate", "difficult", "very_difficult")) {
  if (diff_chr %in% names(silencing_tbl)) {
    writeLines(text = paste("  ", diff_chr, ":",
                           silencing_tbl[diff_chr], "guides",
                           paste0("(", round(100 * silencing_tbl[diff_chr] / nrow(guides_df), 1), "%)")),
               con = summary_file_con)
  }
}

writeLines(text = "\nOverall quality tiers:", con = summary_file_con)
writeLines(text = paste("  Perfect (3 flags + easy silencing):",
                       sum(guides_df$quality_flags_passed_int == 3 & 
                           guides_df$silencing_difficulty_score_int == 1),
                       "guides"),
           con = summary_file_con)
writeLines(text = paste("  High (2-3 flags + easy/moderate silencing):",
                       sum(guides_df$quality_flags_passed_int >= 2 & 
                           guides_df$silencing_difficulty_score_int <= 2),
                       "guides"),
           con = summary_file_con)
writeLines(text = paste("  Moderate (1-2 flags):",
                       sum(guides_df$quality_flags_passed_int >= 1 & 
                           guides_df$quality_flags_passed_int < 3),
                       "guides"),
           con = summary_file_con)
writeLines(text = paste("  Low (0 flags or difficult silencing):",
                       sum(guides_df$quality_flags_passed_int == 0 | 
                           guides_df$silencing_difficulty_score_int >= 3),
                       "guides"),
           con = summary_file_con)

# ============================================================================
# GC content analysis
# ============================================================================
write_section("GC CONTENT ANALYSIS")

writeLines(text = paste("Optimal range:", GC_MIN_PERCENT_dbl, "-", 
                       GC_MAX_PERCENT_dbl, "%"),
           con = summary_file_con)
writeLines(text = paste("Guides in optimal range:",
                       sum(guides_df$gc_in_optimal_range_lgl),
                       paste0("(", round(100 * mean(guides_df$gc_in_optimal_range_lgl), 1), "%)")),
           con = summary_file_con)
writeLines(text = paste("Guides below", GC_MIN_PERCENT_dbl, "%:",
                       sum(guides_df$gc_percent_dbl < GC_MIN_PERCENT_dbl),
                       paste0("(", round(100 * mean(guides_df$gc_percent_dbl < GC_MIN_PERCENT_dbl), 1), "%)")),
           con = summary_file_con)
writeLines(text = paste("Guides above", GC_MAX_PERCENT_dbl, "%:",
                       sum(guides_df$gc_percent_dbl > GC_MAX_PERCENT_dbl),
                       paste0("(", round(100 * mean(guides_df$gc_percent_dbl > GC_MAX_PERCENT_dbl), 1), "%)")),
           con = summary_file_con)

writeLines(text = "\nGC% distribution:", con = summary_file_con)
gc_summary <- summary(guides_df$gc_percent_dbl)
for (i in seq_along(gc_summary)) {
  writeLines(text = paste("  ", names(gc_summary)[i], ":",
                         round(gc_summary[i], 1), "%"),
             con = summary_file_con)
}

# ============================================================================
# Homopolymer analysis
# ============================================================================
write_section("HOMOPOLYMER ANALYSIS")

writeLines(text = paste("Poly-T threshold:", POLYT_THRESHOLD_bp_int, "bp"),
           con = summary_file_con)
writeLines(text = paste("Guides with poly-T problems:",
                       sum(guides_df$has_polyT_problem_lgl),
                       paste0("(", round(100 * mean(guides_df$has_polyT_problem_lgl), 1), "%)")),
           con = summary_file_con)

writeLines(text = "\nGeneral homopolymer threshold:", con = summary_file_con)
writeLines(text = paste("  ", HOMOPOLYMER_THRESHOLD_bp_int, "bp"),
           con = summary_file_con)
writeLines(text = paste("Guides with homopolymer problems:",
                       sum(guides_df$has_homopolymer_problem_lgl),
                       paste0("(", round(100 * mean(guides_df$has_homopolymer_problem_lgl), 1), "%)")),
           con = summary_file_con)

writeLines(text = "\nHomopolymer nucleotide distribution:", 
           con = summary_file_con)
homopoly_nuc_tbl <- table(guides_df$homopolymer_nucleotide_chr, useNA = "ifany")
for (nuc_chr in c("A", "T", "G", "C")) {
  if (nuc_chr %in% names(homopoly_nuc_tbl)) {
    writeLines(text = paste("  ", nuc_chr, ":",
                           homopoly_nuc_tbl[nuc_chr], "guides",
                           paste0("(", round(100 * homopoly_nuc_tbl[nuc_chr] / sum(!is.na(guides_df$homopolymer_nucleotide_chr)), 1), "%)")),
               con = summary_file_con)
  }
}

# ============================================================================
# Per-gene summary
# ============================================================================
write_section("PER-GENE SUMMARY")

writeLines(text = "Top guide per gene:\n", con = summary_file_con)

# Get top guide for each gene
top1_per_gene_df <- guides_df[guides_df$gene_rank_int == 1, ]
top1_per_gene_df <- top1_per_gene_df[order(top1_per_gene_df$gene_name), ]

for (i in seq_len(nrow(top1_per_gene_df))) {
  gene_chr <- top1_per_gene_df$gene_name[i]
  
  writeLines(text = paste("Gene:", gene_chr), con = summary_file_con)
  writeLines(text = paste("  Sequence:", top1_per_gene_df$guide_sequence[i]),
             con = summary_file_con)
  writeLines(text = paste("  Quality score:",
                         round(top1_per_gene_df$overall_quality_score_dbl[i], 2)),
             con = summary_file_con)
  writeLines(text = paste("  Quality flags passed:",
                         top1_per_gene_df$quality_flags_passed_int[i], "/ 3"),
             con = summary_file_con)
  writeLines(text = paste("  GC%:", round(top1_per_gene_df$gc_percent_dbl[i], 1)),
             con = summary_file_con)
  writeLines(text = paste("  Silencing difficulty:",
                         top1_per_gene_df$silencing_difficulty[i]),
             con = summary_file_con)
  writeLines(text = paste("  Total guides available for gene:",
                         guides_per_gene_tbl[gene_chr]),
             con = summary_file_con)
  writeLines(text = "", con = summary_file_con)
}

# ============================================================================
# Recommendations
# ============================================================================
write_section("RECOMMENDATIONS")

writeLines(text = paste("Total high-quality guides available:",
                       sum(guides_df$quality_flags_passed_int >= 2 &
                           guides_df$silencing_difficulty_score_int <= 2)),
           con = summary_file_con)

writeLines(text = "\nGenes with excellent options (>= 10 perfect/high quality guides):",
           con = summary_file_con)

for (gene_chr in names(guides_per_gene_tbl)) {
  gene_guides_df <- guides_df[guides_df$gene_name == gene_chr, ]
  high_quality_count <- sum(gene_guides_df$quality_flags_passed_int >= 2 &
                            gene_guides_df$silencing_difficulty_score_int <= 2)
  
  if (high_quality_count >= 10) {
    writeLines(text = paste("  ", gene_chr, ":", high_quality_count, "high-quality guides"),
               con = summary_file_con)
  }
}

writeLines(text = "\nGenes with limited options (<= 5 high quality guides):",
           con = summary_file_con)

limited_genes_found <- FALSE
for (gene_chr in names(guides_per_gene_tbl)) {
  gene_guides_df <- guides_df[guides_df$gene_name == gene_chr, ]
  high_quality_count <- sum(gene_guides_df$quality_flags_passed_int >= 2 &
                            gene_guides_df$silencing_difficulty_score_int <= 2)
  
  if (high_quality_count <= 5) {
    writeLines(text = paste("  ", gene_chr, ":", high_quality_count, "high-quality guides"),
               con = summary_file_con)
    limited_genes_found <- TRUE
  }
}

if (!limited_genes_found) {
  writeLines(text = "  None - all genes have good guide availability",
             con = summary_file_con)
}

# ============================================================================
# Footer
# ============================================================================
write_section("OUTPUT FILES")

writeLines(text = paste("Full dataset:", OUTPUT_FULL_PATH_chr),
           con = summary_file_con)
writeLines(text = paste("Top", TOP_N_GUIDES_int, "per gene:", OUTPUT_TOP5_PATH_chr),
           con = summary_file_con)
writeLines(text = paste("This summary:", OUTPUT_SUMMARY_PATH_chr),
           con = summary_file_con)

# Close file connection
close(summary_file_con)

cat("=== Output 3: Summary Statistics ===\n")
cat("Written to:", OUTPUT_SUMMARY_PATH_chr, "\n")
cat("\n")

cat("=== ALL OUTPUTS COMPLETE ===\n")
cat("1. Full dataset with quality metrics:", OUTPUT_FULL_PATH_chr, "\n")
cat("2. Top", TOP_N_GUIDES_int, "guides per gene:", OUTPUT_TOP5_PATH_chr, "\n")
cat("3. Summary statistics report:", OUTPUT_SUMMARY_PATH_chr, "\n")
cat("\n")
