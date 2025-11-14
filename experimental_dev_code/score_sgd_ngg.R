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
