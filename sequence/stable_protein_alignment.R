#!/usr/bin/env Rscript
# stable_msa_conservation.R
# Align stable protein sequences and calculate conservation relative to S. cerevisiae

library(Biostrings)
library(DECIPHER)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

INPUT_DIR <- "~/data/protein_files"
OUTPUT_DIR <- "~/data/protein_files/alignments"

GENE_NAMES <- c("ORC1", "ORC2", "ORC3", "ORC4", "ORC5", "ORC6", "TOA2")

INPUT_PREFIX <- "stable_"
OUTPUT_PREFIX <- "stable_"

# Reference organism for conservation (must match FASTA header pattern)
REFERENCE_PATTERN <- "^Scer_"

# ==============================================================================
# SETUP
# ==============================================================================

if (!dir.exists(OUTPUT_DIR)) {
    dir.create(path = OUTPUT_DIR, recursive = TRUE)
}

# Build file paths
input_files <- file.path(INPUT_DIR, paste0(INPUT_PREFIX, GENE_NAMES, ".fasta"))
output_aligned <- file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, GENE_NAMES, "_aligned.fasta"))
output_conservation <- file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, GENE_NAMES, "_conservation.tsv"))

# Check inputs exist
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
    stop("Missing input files:\n", paste(missing_files, collapse = "\n"))
}

# Pre-allocate summary table
summary_df <- data.frame(
    gene = GENE_NAMES,
    n_sequences = integer(length(GENE_NAMES)),
    alignment_length = integer(length(GENE_NAMES)),
    scer_length = integer(length(GENE_NAMES)),
    mean_conservation = numeric(length(GENE_NAMES)),
    stringsAsFactors = FALSE
)

cat("=== STABLE MSA & CONSERVATION ===\n")
cat("Genes:", length(GENE_NAMES), "\n")
cat("Input:", INPUT_DIR, "\n")
cat("Output:", OUTPUT_DIR, "\n\n")

# ==============================================================================
# PROCESSING - ALIGNMENT AND CONSERVATION
# ==============================================================================

for (i in seq_along(GENE_NAMES)) {

    gene <- GENE_NAMES[i]
    cat("[", i, "/", length(GENE_NAMES), "] ", gene, "\n", sep = "")

    # --- Load sequences ---
    seqs <- Biostrings::readAAStringSet(filepath = input_files[i])
    n_seqs <- length(seqs)
    cat("  Loaded:", n_seqs, "sequences\n")

    # --- Align ---
    cat("  Aligning...")
    aligned <- DECIPHER::AlignSeqs(myXStringSet = seqs, verbose = FALSE)
    cat(" done\n")

    # Save alignment
    Biostrings::writeXStringSet(x = aligned, filepath = output_aligned[i])

    aln_length <- unique(Biostrings::width(aligned))
    cat("  Alignment:", n_seqs, "x", aln_length, "\n")

    # --- Find reference sequence ---
    seq_names <- names(aligned)
    ref_idx <- grep(pattern = REFERENCE_PATTERN, x = seq_names)

    if (length(ref_idx) == 0) {
        cat("  WARNING: Reference not found, skipping conservation\n\n")
        summary_df$n_sequences[i] <- n_seqs
        summary_df$alignment_length[i] <- aln_length
        summary_df$scer_length[i] <- NA
        summary_df$mean_conservation[i] <- NA
        next
    }

    if (length(ref_idx) > 1) {
        cat("  WARNING: Multiple references, using first\n")
        ref_idx <- ref_idx[1]
    }

    # --- Conservation calculation ---
    # Convert to matrix (rows = sequences, cols = positions)
    aln_matrix <- as.matrix(aligned)
    ref_seq <- aln_matrix[ref_idx, ]

    # Identify ungapped reference positions
    ungapped_cols <- which(ref_seq != "-")
    scer_length <- length(ungapped_cols)

    # Pre-allocate conservation profile
    conservation_df <- data.frame(
        position = seq_len(scer_length),
        residue = ref_seq[ungapped_cols],
        conservation_pct = numeric(scer_length),
        stringsAsFactors = FALSE
    )

    # Calculate conservation at each reference position
    # Exclude reference row from calculation
    other_rows <- seq_len(n_seqs)[-ref_idx]
    other_matrix <- aln_matrix[other_rows, , drop = FALSE]
    n_other <- length(other_rows)

    for (j in seq_along(ungapped_cols)) {
        col_idx <- ungapped_cols[j]
        ref_residue <- ref_seq[col_idx]
        col_values <- other_matrix[, col_idx]

        # Count non-gap sequences
        non_gap <- col_values != "-"
        n_non_gap <- sum(non_gap)

        if (n_non_gap > 0) {
            n_identical <- sum(col_values[non_gap] == ref_residue)
            conservation_df$conservation_pct[j] <- (n_identical / n_non_gap) * 100
        } else {
            conservation_df$conservation_pct[j] <- NA
        }
    }

    # Save conservation profile
    utils::write.table(
        x = conservation_df,
        file = output_conservation[i],
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    mean_cons <- mean(conservation_df$conservation_pct, na.rm = TRUE)
    cat("  Conservation: mean =", round(mean_cons, 1), "% across", scer_length, "positions\n")

    # Update summary
    summary_df$n_sequences[i] <- n_seqs
    summary_df$alignment_length[i] <- aln_length
    summary_df$scer_length[i] <- scer_length
    summary_df$mean_conservation[i] <- mean_cons

    cat("\n")
}

# ==============================================================================
# OUTPUT - SUMMARY
# ==============================================================================

summary_path <- file.path(OUTPUT_DIR, paste0(OUTPUT_PREFIX, "summary.tsv"))
utils::write.table(
    x = summary_df,
    file = summary_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

cat("=== COMPLETE ===\n")
cat("Summary saved to:", summary_path, "\n\n")
print(summary_df)
