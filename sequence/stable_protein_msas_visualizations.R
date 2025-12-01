#!/usr/bin/env Rscript
# stable_msa_visualization.R
# Generate zoomed MSA visualizations at specified S. cerevisiae positions
# Requires devtools-mediated installation of ggtree then ggmsa.

library(Biostrings)
library(ggmsa)
library(ggplot2)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

ALIGNMENT_DIR <- "~/data/protein_files/alignments"
PLOT_OUTPUT_DIR <- "~/data/protein_files/alignments/plots"

INPUT_PREFIX <- "stable_"

GENE_NAMES <- c("ORC1", "ORC2", "ORC3", "ORC4", "ORC5", "ORC6", "TOA2")

# Reference organism pattern (must match FASTA header)
REFERENCE_PATTERN <- "^Scer_"

# ==============================================================================
# REGIONS OF INTEREST
# ==============================================================================
# Positions in UNGAPPED S. cerevisiae coordinates
# Script maps these to gapped alignment positions automatically

REGIONS_OF_INTEREST <- list(
    ORC1 = c(495),
    ORC2 = c(),
    ORC3 = c(481),
    ORC4 = c(225),
    ORC5 = c(104),
    ORC6 = c(305),
    TOA2 = c(16)
)

# Window: residues on each side of center position
ZOOM_WINDOW <- 15

# ==============================================================================
# PLOT SETTINGS
# ==============================================================================

COLOR_SCHEME <- "Chemistry_AA"
CHAR_WIDTH <- 0.5
SHOW_SEQ_NAMES <- TRUE
SHOW_LOGO <- FALSE

PLOT_WIDTH <- 12
PLOT_HEIGHT <- 8

SAVE_PNG <- TRUE
SAVE_SVG <- TRUE

# ==============================================================================
# SETUP
# ==============================================================================

if (!dir.exists(PLOT_OUTPUT_DIR)) {
    dir.create(path = PLOT_OUTPUT_DIR, recursive = TRUE)
}

cat("=== STABLE MSA VISUALIZATION ===\n")
cat("Input:", ALIGNMENT_DIR, "\n")
cat("Output:", PLOT_OUTPUT_DIR, "\n\n")

# Count total plots to generate
total_positions <- sum(sapply(REGIONS_OF_INTEREST, length))
cat("Positions to visualize:", total_positions, "\n\n")

# Track results
results <- list()

# ==============================================================================
# PROCESSING
# ==============================================================================

for (gene in GENE_NAMES) {

    positions <- REGIONS_OF_INTEREST[[gene]]

    if (length(positions) == 0) {
        cat("Skipping", gene, "(no positions)\n")
        next
    }

    cat("[", gene, "] Positions:", paste(positions, collapse = ", "), "\n", sep = "")

    # Load alignment
    alignment_path <- file.path(ALIGNMENT_DIR, paste0(INPUT_PREFIX, gene, "_aligned.fasta"))

    if (!file.exists(alignment_path)) {
        cat("  WARNING: Alignment not found, skipping\n\n")
        for (pos in positions) {
            results[[length(results) + 1]] <- list(
                gene = gene, position = pos, status = "alignment_missing"
            )
        }
        next
    }

    aligned <- Biostrings::readAAStringSet(filepath = alignment_path)

    # Find reference sequence
    seq_names <- names(aligned)
    ref_idx <- grep(pattern = REFERENCE_PATTERN, x = seq_names)

    if (length(ref_idx) == 0) {
        cat("  WARNING: Reference not found, skipping\n\n")
        for (pos in positions) {
            results[[length(results) + 1]] <- list(
                gene = gene, position = pos, status = "reference_missing"
            )
        }
        next
    }

    if (length(ref_idx) > 1) {
        cat("  WARNING: Multiple references, using first\n")
        ref_idx <- ref_idx[1]
    }

    # Extract reference sequence for position mapping
    ref_seq <- strsplit(as.character(aligned[ref_idx]), "")[[1]]
    aln_length <- length(ref_seq)

    # Build ungapped-to-gapped position map (once per gene)
    ungapped_pos <- 0
    pos_map <- integer(length(ref_seq))  # ungapped -> gapped

    for (i in seq_along(ref_seq)) {
        if (ref_seq[i] != "-") {
            ungapped_pos <- ungapped_pos + 1
            pos_map[ungapped_pos] <- i
        }
    }
    ref_length <- ungapped_pos  # total ungapped positions

    # Process each position
    for (pos in positions) {

        cat("  Position", pos, "... ")

        # Check bounds
        if (pos < 1 || pos > ref_length) {
            cat("out of range (max:", ref_length, ")\n")
            results[[length(results) + 1]] <- list(
                gene = gene, position = pos, status = "out_of_range"
            )
            next
        }

        # Map to alignment position
        aln_pos <- pos_map[pos]
        residue <- ref_seq[aln_pos]

        # Calculate window
        start_pos <- max(1, aln_pos - ZOOM_WINDOW)
        end_pos <- min(aln_length, aln_pos + ZOOM_WINDOW)

        cat("aln:", start_pos, "-", end_pos, "| residue:", residue, "\n")

        # Generate plot
        p <- ggmsa::ggmsa(
            msa = aligned,
            start = start_pos,
            end = end_pos,
            color = COLOR_SCHEME,
            char_width = CHAR_WIDTH,
            seq_name = SHOW_SEQ_NAMES,
            consensus_views = SHOW_LOGO,
            disagreement = FALSE,
            use_dot = FALSE,
            by_conservation = FALSE,
            none_bg = FALSE
        ) +
            ggplot2::ggtitle(paste0(gene, " - Position ", pos, " (", residue, ")")) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold")
            )

        # Save
        output_base <- paste0(INPUT_PREFIX, gene, "_pos", pos)

        if (SAVE_PNG) {
            png_path <- file.path(PLOT_OUTPUT_DIR, paste0(output_base, ".png"))
            ggplot2::ggsave(
                filename = png_path,
                plot = p,
                width = PLOT_WIDTH,
                height = PLOT_HEIGHT,
                dpi = 300
            )
        }

        if (SAVE_SVG) {
            svg_path <- file.path(PLOT_OUTPUT_DIR, paste0(output_base, ".svg"))
            ggplot2::ggsave(
                filename = svg_path,
                plot = p,
                width = PLOT_WIDTH,
                height = PLOT_HEIGHT
            )
        }

        results[[length(results) + 1]] <- list(
            gene = gene,
            position = pos,
            status = "complete",
            file = output_base
        )
    }

    cat("\n")
}

# ==============================================================================
# OUTPUT - SUMMARY
# ==============================================================================

summary_df <- do.call(rbind, lapply(results, as.data.frame))

summary_path <- file.path(PLOT_OUTPUT_DIR, paste0(INPUT_PREFIX, "visualization_summary.tsv"))
utils::write.table(
    x = summary_df,
    file = summary_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

cat("=== COMPLETE ===\n")
cat("Plots generated:", sum(summary_df$status == "complete"), "/", nrow(summary_df), "\n")
cat("Summary:", summary_path, "\n\n")
print(summary_df)
