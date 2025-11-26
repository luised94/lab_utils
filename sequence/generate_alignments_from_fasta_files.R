# ==============================================================================
# MSA Generation and Conservation Aanlysis.
# ==============================================================================
# Purpose: Generate multiple sequence alignments for ORC complex and TFIIA
#          proteins, calculate conservation statistics, and create zoomed
#          visualizations at specified residue positions
#
# Input:   Filtered UniRef50 FASTA files and model organism sequences
# Output:  Aligned FASTA files, conservation profiles, zoomed region plots
# Date: 2025-11-25
# ==============================================================================

# === PACKAGE LOADING ==========================================================

library(Biostrings)  # Bioconductor - FASTA I/O, sequence manipulation
library(DECIPHER)    # Bioconductor - MUSCLE alignment (AlignSeqs function)

# === FILE PATHS ===============================================================

# Input directories
INPUT_DIR_path <- "~/data/protein_files"
ALIGNMENT_DIR_path <- "~/data/protein_files/alignments"

# Output directory (all outputs saved here, flat structure)
OUTPUT_DIR_path <- "~/data/protein_files/alignments"

# === GENES TO PROCESS =========================================================

GENE_NAMES_chr <- c("ORC1", "ORC2", "ORC3", "ORC4", "ORC5", "ORC6", 
                    "TOA1", "TOA2")

# === ALIGNMENT CONFIGURATION ==================================================

# Alignment method (ClustalOmega via msa package)
ALIGNMENT_METHOD_chr <- "MUSCLE"

# Preserve input sequence order in alignments
ALIGNMENT_ORDER_chr <- "input"

# Force re-alignment even if cached files exist?
FORCE_REALIGNMENT_lgl <- FALSE

# Cache aligned FASTA files for faster re-runs
USE_ALIGNMENT_CACHE_lgl <- TRUE

# === REGIONS OF INTEREST ======================================================
# Positions refer to UNGAPPED S. cerevisiae reference sequence coordinates
# The script will map these to gapped alignment positions automatically
# 
# UPDATE THESE with your actual positions of interest
# Leave as c() to skip visualization for that gene

REGIONS_OF_INTEREST_lst <- list(
  ORC1 = c(495),   # Example positions - UPDATE AS NEEDED
  ORC2 = c(),         # Example positions - UPDATE AS NEEDED
  ORC3 = c(481),         # Example positions - UPDATE AS NEEDED
  ORC4 = c(225),              # Example positions - UPDATE AS NEEDED
  ORC5 = c(104),         # Example positions - UPDATE AS NEEDED
  ORC6 = c(305),                 # No visualizations for this gene
  TOA1 = c(49),         # Example positions - UPDATE AS NEEDED
  TOA2 = c(16)           # Example positions - UPDATE AS NEEDED
)

# Window size for zoomed plots:  residues around center position
ZOOM_WINDOW_RESIDUES_int <- 10

# === CONSERVATION ANALYSIS ====================================================

# Calculate conservation as percent identity per position
# Profile will include all S. cerevisiae positions (ungapped coordinates)
CONSERVATION_METHOD_chr <- "percent_identity"

# Generate conservation profiles for both datasets
CALCULATE_MODEL_ORGS_CONSERVATION_lgl <- TRUE
CALCULATE_UNIREF50_CONSERVATION_lgl <- TRUE

# === UNIREF50 SUBSETTING FOR VISUALIZATION ===================================

# Number of sequences to show (including S. cerevisiae)
UNIREF50_PLOT_N_SEQUENCES_int <- 15

# Subset UniRef50 alignments to top N sequences for cleaner plots
# Ranking is by percent identity to S. cerevisiae
SUBSET_UNIREF50_FOR_PLOTS_lgl <- TRUE

# Always include S. cerevisiae as first sequence
UNIREF50_KEEP_REFERENCE_lgl <- TRUE

# === S. CEREVISIAE REFERENCE IDENTIFICATION ==================================

# Pattern to identify S. cerevisiae sequences in alignments
# Model organisms format: "Scer_GENE|..." 
# UniRef50 format: "Sac Cer GENE|..."
SCER_PATTERN_MODEL_ORGS_chr <- "^Scer_"
SCER_PATTERN_UNIREF50_chr <- "^Sac Cer"

# ==============================================================================
# End of Configuration
# ==============================================================================

# ==============================================================================
# CHUNK 1: ALIGNMENT GENERATION
# ==============================================================================
# Generate multiple sequence alignments using DECIPHER::AlignSeqs (MUSCLE)
# Process both model organism and UniRef50 datasets for all genes

cat("\n=== STARTING ALIGNMENT GENERATION ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Initialize tracking
alignment_summary_df <- data.frame(
  gene = character(),
  dataset = character(),
  n_sequences = integer(),
  alignment_length = integer(),
  status = character(),
  stringsAsFactors = FALSE
)

# Process each gene
for (gene_chr in GENE_NAMES_chr) {

  cat("Processing gene:", gene_chr, "\n")

  # === MODEL ORGANISMS ALIGNMENT ===
  cat("  [1/2] Model organisms...\n")

  # Construct file paths
  input_fasta_path <- file.path(INPUT_DIR_path, 
                                 paste0(gene_chr, "_model_organisms.fasta"))
  output_fasta_path <- file.path(OUTPUT_DIR_path,
                                  paste0(gene_chr, "_model_orgs_aligned.fasta"))

  # Check if input exists
  if (!file.exists(input_fasta_path)) {
    cat("    WARNING: Input file not found, skipping\n")
    alignment_summary_df <- rbind(alignment_summary_df,
                                  data.frame(gene = gene_chr,
                                           dataset = "model_orgs",
                                           n_sequences = NA,
                                           alignment_length = NA,
                                           status = "input_missing"))
  } else {

    # Check cache
    skip_alignment_lgl <- FALSE
    if (USE_ALIGNMENT_CACHE_lgl && !FORCE_REALIGNMENT_lgl && 
        file.exists(output_fasta_path)) {
      cat("    Using cached alignment\n")
      skip_alignment_lgl <- TRUE
    }

    if (!skip_alignment_lgl) {
      # Load sequences
      seqs_AAStringSet <- Biostrings::readAAStringSet(filepath = input_fasta_path)
      cat("    Loaded", length(seqs_AAStringSet), "sequences\n")

      # Align using DECIPHER
      cat("    Running MUSCLE alignment...\n")
      aligned_AAStringSet <- DECIPHER::AlignSeqs(myXStringSet = seqs_AAStringSet,
                                                  verbose = FALSE)

      # Save aligned sequences
      Biostrings::writeXStringSet(x = aligned_AAStringSet,
                                  filepath = output_fasta_path)
      cat("    Saved aligned FASTA\n")
    } else {
      # Load cached alignment for reporting
      aligned_AAStringSet <- Biostrings::readAAStringSet(filepath = output_fasta_path)
    }

    # Report dimensions
    n_seqs_int <- length(aligned_AAStringSet)
    aln_length_int <- unique(width(aligned_AAStringSet))
    cat("    Alignment dimensions:", n_seqs_int, "sequences x", 
        aln_length_int, "positions\n")

    # Record summary
    alignment_summary_df <- rbind(alignment_summary_df,
                                  data.frame(gene = gene_chr,
                                           dataset = "model_orgs",
                                           n_sequences = n_seqs_int,
                                           alignment_length = aln_length_int,
                                           status = "complete"))
  }

  # === UNIREF50 ALIGNMENT ===
  cat("  [2/2] UniRef50...\n")

  # Construct file paths
  input_fasta_path <- file.path(ALIGNMENT_DIR_path,
                                 paste0(gene_chr, "_uniref50_filtered.fasta"))
  output_fasta_path <- file.path(OUTPUT_DIR_path,
                                  paste0(gene_chr, "_uniref50_aligned.fasta"))

  # Check if input exists
  if (!file.exists(input_fasta_path)) {
    cat("    WARNING: Input file not found, skipping\n")
    alignment_summary_df <- rbind(alignment_summary_df,
                                  data.frame(gene = gene_chr,
                                           dataset = "uniref50",
                                           n_sequences = NA,
                                           alignment_length = NA,
                                           status = "input_missing"))
  } else {

    # Check cache
    skip_alignment_lgl <- FALSE
    if (USE_ALIGNMENT_CACHE_lgl && !FORCE_REALIGNMENT_lgl && 
        file.exists(output_fasta_path)) {
      cat("    Using cached alignment\n")
      skip_alignment_lgl <- TRUE
    }

    if (!skip_alignment_lgl) {
      # Load sequences
      seqs_AAStringSet <- Biostrings::readAAStringSet(filepath = input_fasta_path)
      cat("    Loaded", length(seqs_AAStringSet), "sequences\n")

      # Align using DECIPHER
      cat("    Running MUSCLE alignment...\n")
      aligned_AAStringSet <- DECIPHER::AlignSeqs(myXStringSet = seqs_AAStringSet,
                                                  verbose = FALSE)

      # Save aligned sequences
      Biostrings::writeXStringSet(x = aligned_AAStringSet,
                                  filepath = output_fasta_path)
      cat("    Saved aligned FASTA\n")
    } else {
      # Load cached alignment for reporting
      aligned_AAStringSet <- Biostrings::readAAStringSet(filepath = output_fasta_path)
    }

    # Report dimensions
    n_seqs_int <- length(aligned_AAStringSet)
    aln_length_int <- unique(width(aligned_AAStringSet))
    cat("    Alignment dimensions:", n_seqs_int, "sequences x", 
        aln_length_int, "positions\n")

    # Record summary
    alignment_summary_df <- rbind(alignment_summary_df,
                                  data.frame(gene = gene_chr,
                                           dataset = "uniref50",
                                           n_sequences = n_seqs_int,
                                           alignment_length = aln_length_int,
                                           status = "complete"))
  }

  cat("\n")
}

# Save alignment summary
summary_path <- file.path(OUTPUT_DIR_path, "alignment_summary.tsv")
write.table(x = alignment_summary_df,
            file = summary_path,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

cat("=== ALIGNMENT GENERATION COMPLETE ===\n")
cat("Summary saved to:", summary_path, "\n")
print(alignment_summary_df)
cat("\n")
