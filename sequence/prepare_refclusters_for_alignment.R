# prepare_refclusters_for_alignment.R
# Purpose:
#   Filter UniRef50 sequences by length, duplicates, and species representation
#   to prepare for msas.
# Date: 2025-11-25

# === REQUIRED PACKAGES ===
library(Biostrings)

# === PATHS ===
INPUT_DIR_path <- "~/data/protein_files"
OUTPUT_DIR_path <- "~/data/protein_files/alignments"

# === GENES TO PROCESS ===
GENE_NAMES_chr <- c("ORC1", "ORC2", "ORC3", "ORC4", "ORC5", "ORC6", "TOA1", "TOA2")

# === FILTERING CRITERIA ===
MIN_SEQUENCE_LENGTH_int <- 50  # minimum amino acids
KEEP_ONE_PER_SPECIES_lgl <- TRUE

# === REGIONS OF INTEREST (Manual specification) ===
# Positions refer to S. cerevisiae reference sequence (ungapped)
# TO BE UPDATED after alignment inspection

REGIONS_OF_INTEREST_lst <- list(
  ORC1 = c(100, 250, 500),  # Example - UPDATE AFTER INSPECTION
  ORC2 = c(150, 300),        # Example - UPDATE AFTER INSPECTION
  ORC3 = c(200, 400),        # Example - UPDATE AFTER INSPECTION
  ORC4 = c(180),             # Example - UPDATE AFTER INSPECTION
  ORC5 = c(220, 450),        # Example - UPDATE AFTER INSPECTION
  ORC6 = c(130, 280),        # Example - UPDATE AFTER INSPECTION
  TOA1 = c(100, 200),        # Example - UPDATE AFTER INSPECTION
  TOA2 = c(90, 180)          # Example - UPDATE AFTER INSPECTION
)

# Window around each position for zoomed plots
ZOOM_WINDOW_RESIDUES_int <- 10  # +/- residues around center position

# === ALIGNMENT CONFIGURATION ===
ALIGNMENT_METHOD_chr <- "ClustalOmega"
ALIGNMENT_ORDER_chr <- "input"
FORCE_REALIGNMENT_lgl <- FALSE

# === VISUALIZATION CONFIGURATION ===
PLOT_COLOR_SCHEME_chr <- "Chemistry_AA"
PLOT_FONT_chr <- "TimesNewRoman"
SHOW_SEQUENCE_LOGO_lgl <- TRUE
PLOT_WIDTH_inches <- 12
PLOT_HEIGHT_inches <- 8

# Create output directory if needed
if (!dir.exists(OUTPUT_DIR_path)) {
  dir.create(OUTPUT_DIR_path, recursive = TRUE)
  cat("Created output directory:", OUTPUT_DIR_path, "\n")
}

# === CHUNK 1.1: LOAD UNIREF50 FASTA FILES ===
cat("\n=== Loading UniRef50 FASTA files ===\n")

# Initialize storage list
sequences_raw_lst <- list()

# Load each gene's UniRef50 sequences
for (gene_chr in GENE_NAMES_chr) {

  # Construct file path
  fasta_file_path <- file.path(INPUT_DIR_path, paste0(gene_chr, "_uniref50.fasta"))

  # Check file exists
  if (!file.exists(fasta_file_path)) {
    cat("WARNING: File not found:", fasta_file_path, "\n")
    next
  }

  # Read FASTA file
  seqs_AAStringSet <- readAAStringSet(filepath = fasta_file_path)

  # Extract information
  seq_ids_chr <- names(seqs_AAStringSet)
  sequences_chr <- as.character(seqs_AAStringSet)
  lengths_int <- width(seqs_AAStringSet)

  # Create data frame
  sequences_raw_lst[[gene_chr]] <- data.frame(
    gene = gene_chr,
    seq_id = seq_ids_chr,
    header_full = seq_ids_chr,  # Keep full header for parsing
    sequence = sequences_chr,
    length = lengths_int,
    stringsAsFactors = FALSE
  )

  cat(sprintf("  %s: Loaded %d sequences (length range: %d-%d aa)\n",
              gene_chr,
              length(seqs_AAStringSet),
              min(lengths_int),
              max(lengths_int)))
}

# Verify all genes loaded
cat("\n=== Loading Summary ===\n")
cat("Genes successfully loaded:", length(sequences_raw_lst), "of", length(GENE_NAMES_chr), "\n")

# Print sample headers to inspect format (first 3 from ORC1)
if ("ORC1" %in% names(sequences_raw_lst)) {
  cat("\n=== Sample Headers (ORC1, first 3) ===\n")
  sample_headers_chr <- head(sequences_raw_lst[["ORC1"]]$header_full, 3)
  for (i in seq_along(sample_headers_chr)) {
    cat(sprintf("[%d] %s\n", i, sample_headers_chr[i]))
  }
}
