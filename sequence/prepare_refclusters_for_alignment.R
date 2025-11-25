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


# === CHUNK 1.2: PARSE HEADERS TO EXTRACT METADATA ===
cat("\n=== Parsing sequence headers ===\n")

# Parse headers for all genes
for (gene_chr in names(sequences_raw_lst)) {

  cat(sprintf("  Parsing %s headers...\n", gene_chr))

  df <- sequences_raw_lst[[gene_chr]]
  n_seqs <- nrow(df)

  # Initialize new columns
  df$accession <- character(n_seqs)
  df$review_status <- character(n_seqs)
  df$organism_full <- character(n_seqs)
  df$organism_short <- character(n_seqs)
  df$taxid <- character(n_seqs)

  # Parse each header
  for (i in 1:n_seqs) {
     # Initialize result with NA values
     result_lst <- list(
       accession = NA_character_,
       review_status = NA_character_,
       organism_full = NA_character_,
       organism_short = NA_character_,
       taxid = NA_character_
     )

    # Get the header information.
    header_chr <- df$header_full[i]

     # Extract accession (second field in pipe-delimited start)
     # Format: sp|P54784|ORC1_YEAST or tr|J6EQC4|J6EQC4_SACK1
     accession_match <- regexpr("^[a-z]+\\|([A-Z0-9]+)\\|", header_chr, perl = TRUE)
     if (accession_match[1] > 0) {
       accession_full <- regmatches(header_chr, accession_match)
       result_lst$accession <- sub("^[a-z]+\\|([A-Z0-9]+)\\|.*", "\\1", accession_full)
     }

     # Extract review status (sp or tr)
     review_match <- regexpr("^(sp|tr)\\|", header_chr, perl = TRUE)
     if (review_match[1] > 0) {
       result_lst$review_status <- sub("^(sp|tr)\\|.*", "\\1", header_chr)
     }

     # Extract full organism (between OS= and OX=)
     os_match <- regexpr("OS=([^=]+)OX=", header_chr, perl = TRUE)
     if (os_match[1] > 0) {
       os_full <- regmatches(header_chr, os_match)
       # Clean up: remove OS= and OX=, trim whitespace
       organism_full <- sub("OS=(.+)OX=", "\\1", os_full)
       organism_full <- gsub("^\\s+|\\s+$", "", organism_full)  # trim
       result_lst$organism_full <- organism_full

       # Extract short name (first two words, typically genus + species)
       # Stop at parenthesis if present
       organism_short <- sub("\\s*\\(.*", "", organism_full)  # remove parenthetical
       words <- strsplit(organism_short, "\\s+")[[1]]
       result_lst$organism_short <- paste(words[1:min(2, length(words))], collapse = " ")
     }

     # Extract taxonomy ID (digits after OX=)
     ox_match <- regexpr("OX=([0-9]+)", header_chr, perl = TRUE)
     if (ox_match[1] > 0) {
       ox_full <- regmatches(header_chr, ox_match)
       result_lst$taxid <- sub("OX=([0-9]+).*", "\\1", ox_full)
     }

    df$accession[i] <- result_lst$accession
    df$review_status[i] <- result_lst$review_status
    df$organism_full[i] <- result_lst$organism_full
    df$organism_short[i] <- result_lst$organism_short
    df$taxid[i] <- result_lst$taxid
  }

  # Update the list
  sequences_raw_lst[[gene_chr]] <- df

  # Report parsing success
  n_parsed <- sum(!is.na(df$taxid))
  cat(sprintf("    Successfully parsed %d/%d sequences (%.1f%%)\n",
              n_parsed, n_seqs, 100 * n_parsed / n_seqs))
}

# Print summary statistics
cat("\n=== Parsing Summary ===\n")
cat("Review status distribution:\n")
for (gene_chr in names(sequences_raw_lst)) {
  df <- sequences_raw_lst[[gene_chr]]
  n_sp <- sum(df$review_status == "sp", na.rm = TRUE)
  n_tr <- sum(df$review_status == "tr", na.rm = TRUE)
  cat(sprintf("  %s: sp=%d, tr=%d\n", gene_chr, n_sp, n_tr))
}

# Show example parsed data (first 3 rows of ORC1)
cat("\n=== Example Parsed Data (ORC1, first 3) ===\n")
if ("ORC1" %in% names(sequences_raw_lst)) {
  example_df <- sequences_raw_lst[["ORC1"]][1:3, c("accession", "review_status",
                                                     "organism_short", "taxid", "length")]
  print(example_df, row.names = FALSE)
}

# === CHUNK 1.3: FILTER BY SEQUENCE LENGTH ===
cat("\n=== Filtering by sequence length (>= ", MIN_SEQUENCE_LENGTH_int, " aa) ===\n")

# Storage for filtering report
length_filter_report_df <- data.frame(
  gene = character(),
  n_original = integer(),
  n_below_threshold = integer(),
  n_retained = integer(),
  percent_retained = numeric(),
  stringsAsFactors = FALSE
)

# Filter each gene
sequences_length_filtered_lst <- list()

for (gene_chr in names(sequences_raw_lst)) {

  df <- sequences_raw_lst[[gene_chr]]
  n_original <- nrow(df)

  # Apply length filter
  df_filtered <- subset(df, length >= MIN_SEQUENCE_LENGTH_int)
  n_retained <- nrow(df_filtered)
  n_removed <- n_original - n_retained

  # Store filtered data
  sequences_length_filtered_lst[[gene_chr]] <- df_filtered

  # Add to report
  length_filter_report_df <- rbind(
    length_filter_report_df,
    data.frame(
      gene = gene_chr,
      n_original = n_original,
      n_below_threshold = n_removed,
      n_retained = n_retained,
      percent_retained = round(100 * n_retained / n_original, 1),
      stringsAsFactors = FALSE
    )
  )

  # Console output
  if (n_removed > 0) {
    cat(sprintf("  %s: Removed %d short sequences, kept %d (%.1f%%)\n",
                gene_chr, n_removed, n_retained, 100 * n_retained / n_original))
  } else {
    cat(sprintf("  %s: All %d sequences passed length filter\n", gene_chr, n_original))
  }
}

# Print summary table
cat("\n=== Length Filter Summary ===\n")
print(length_filter_report_df, row.names = FALSE)

# Write report to file
report_file_path <- file.path(OUTPUT_DIR_path, "uniref50_length_filter_report.tsv")
write.table(
  length_filter_report_df,
  file = report_file_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
cat("\nReport saved:", report_file_path, "\n")

# Verify S. cerevisiae retained in all genes
cat("\n=== Verification: S. cerevisiae sequences retained ===\n")
for (gene_chr in names(sequences_length_filtered_lst)) {
  df <- sequences_length_filtered_lst[[gene_chr]]
  scer_present <- any(df$taxid == "559292", na.rm = TRUE)
  cat(sprintf("  %s: %s\n", gene_chr, ifelse(scer_present, "YES", "MISSING!")))
}
