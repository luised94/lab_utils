# ==============================================================================
# MSA Generation and Visualization Pipeline
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
library(DECIPHER)         # Bioconductor - ClustalOmega alignment wrapper
library(ggmsa)       # CRAN/Bioconductor - alignment visualization
library(ggplot2)     # CRAN - plot customization and themes

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
ALIGNMENT_METHOD_chr <- "ClustalOmega"

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

# Subset UniRef50 alignments to top N sequences for cleaner plots
# Ranking is by percent identity to S. cerevisiae
SUBSET_UNIREF50_FOR_PLOTS_lgl <- TRUE

# Number of sequences to show (including S. cerevisiae)
UNIREF50_PLOT_N_SEQUENCES_int <- 15

# Always include S. cerevisiae as first sequence
UNIREF50_KEEP_REFERENCE_lgl <- TRUE

# === VISUALIZATION CONFIGURATION ==============================================

# Color scheme for amino acids
# Options: "Chemistry_AA", "Shapely_AA", "Zappo_AA", "Taylor_AA", "Clustal"
PLOT_COLOR_SCHEME_chr <- "Chemistry_AA"

# Font family for plots
PLOT_FONT_chr <- "sans"  # "sans", "serif", "mono", or specific font names

# Show sequence logo (conservation visualization) above alignment
SHOW_SEQUENCE_LOGO_lgl <- TRUE

# Character width for amino acid letters (smaller = more compact)
PLOT_CHAR_WIDTH_dbl <- 0.5

# Show sequence names on plots
SHOW_SEQUENCE_NAMES_lgl <- TRUE

# Plot dimensions (inches)
PLOT_WIDTH_inches <- 12
PLOT_HEIGHT_inches <- 8

# Output format (always generate PNG, optionally SVG)
SAVE_PNG_lgl <- TRUE
SAVE_SVG_lgl <- FALSE

# === S. CEREVISIAE REFERENCE IDENTIFICATION ==================================

# Pattern to identify S. cerevisiae sequences in alignments
# Model organisms format: "Scer_GENE|..." 
# UniRef50 format: "Sac Cer GENE|..."
SCER_PATTERN_MODEL_ORGS_chr <- "^Scer_"
SCER_PATTERN_UNIREF50_chr <- "^Sac Cer"

# ==============================================================================
# End of Configuration
# ==============================================================================
