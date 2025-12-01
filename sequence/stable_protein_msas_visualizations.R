#!/usr/bin/env Rscript
# stable_msa_visualization.R
# Generate zoomed MSA visualizations at specified S. cerevisiae positions
#
# NOTE: ggmsa produces deprecation warnings about aes_() and "No shared levels"
# warnings related to fill scales. These are internal to ggmsa and can be ignored.
# See: https://github.com/YuLab-SMU/ggmsa/issues
# Requires devtools-mediated installation of ggtree then ggmsa.

library(Biostrings)
library(DECIPHER)
library(ggmsa)
library(ggplot2)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

ALIGNMENT_DIR <- "~/data/protein_files/alignments"
PLOT_OUTPUT_DIR <- "~/data/protein_files/alignments/plots"

INPUT_PREFIX <- "stable_"

GENE_NAMES <- c("ORC1", "ORC2", "ORC3", "ORC4", "ORC5", "ORC6", "TOA2")

# Reference organism pattern (must match FASTA header prefix)
REFERENCE_PATTERN <- "^Scer_"
REFERENCE_ORG <- "Scer"

# ==============================================================================
# SEQUENCE ORDERING CONFIGURATION
# ==============================================================================

# Method: "identity" (% identity to reference) or "hierarchical" (clustering)
ORDERING_METHOD <- "identity"

# Gene used to determine sequence order (applied to all genes)
ORDERING_REFERENCE_GENE <- "ORC1"

# For hierarchical clustering
HCLUST_METHOD <- "average"

# ==============================================================================
# REGIONS OF INTEREST
# ==============================================================================
# Positions in UNGAPPED S. cerevisiae coordinates

REGIONS_OF_INTEREST <- list(
    ORC1 = c(495),
    ORC2 = c(),
    ORC3 = c(481),
    ORC4 = c(225),
    ORC5 = c(104),
    ORC6 = c(305),
    TOA2 = c(16)
)

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
cat("Ordering method:", ORDERING_METHOD, "\n")
cat("Ordering reference gene:", ORDERING_REFERENCE_GENE, "\n\n")

# ==============================================================================
# DETERMINE SEQUENCE ORDER FROM REFERENCE GENE
# ==============================================================================

cat("Computing sequence order from", ORDERING_REFERENCE_GENE, "...\n")

ref_gene_path <- file.path(
    ALIGNMENT_DIR,
    paste0(INPUT_PREFIX, ORDERING_REFERENCE_GENE, "_aligned.fasta")
)

if (!file.exists(ref_gene_path)) {
    stop("Reference gene alignment not found: ", ref_gene_path)
}

ref_aligned <- Biostrings::readAAStringSet(filepath = ref_gene_path)
ref_names <- names(ref_aligned)

# Extract organism codes from headers (format: "Organism_Gene")
ref_organisms <- sub(pattern = "_.*$", replacement = "", x = ref_names)

# Find reference organism index
ref_idx <- grep(pattern = REFERENCE_PATTERN, x = ref_names)
if (length(ref_idx) == 0) {
    stop("Reference organism not found in ", ORDERING_REFERENCE_GENE)
}
if (length(ref_idx) > 1) {
    ref_idx <- ref_idx[1]
}

# Get non-reference indices
other_idx <- setdiff(seq_along(ref_names), ref_idx)
other_organisms <- ref_organisms[other_idx]

# Compute ordering based on method
if (ORDERING_METHOD == "identity") {

    cat("  Computing percent identity to", REFERENCE_ORG, "...\n")

    ref_seq <- as.character(ref_aligned[ref_idx])
    ref_chars <- strsplit(ref_seq, "")[[1]]

    identity_scores <- numeric(length(other_idx))

    for (i in seq_along(other_idx)) {
        idx <- other_idx[i]
        seq_chars <- strsplit(as.character(ref_aligned[idx]), "")[[1]]

        # Compare at positions where both are non-gap
        both_present <- ref_chars != "-" & seq_chars != "-"
        n_compared <- sum(both_present)

        if (n_compared > 0) {
            n_identical <- sum(ref_chars[both_present] == seq_chars[both_present])
            identity_scores[i] <- n_identical / n_compared * 100
        } else {
            identity_scores[i] <- 0
        }
    }

    # Order by identity descending
    order_idx <- order(identity_scores, decreasing = TRUE)
    ordered_organisms <- other_organisms[order_idx]

    cat("  Identity scores:\n")
    for (i in order_idx) {
        cat("    ", other_organisms[i], ": ", round(identity_scores[i], 1), "%\n", sep = "")
    }

} else if (ORDERING_METHOD == "hierarchical") {

    cat("  Computing distance matrix and clustering...\n")

    dist_matrix <- DECIPHER::DistanceMatrix(
        myXStringSet = ref_aligned,
        type = "dist",
        verbose = FALSE
    )

    hc <- stats::hclust(d = as.dist(dist_matrix), method = HCLUST_METHOD)

    # Get leaf order from dendrogram
    dend_order <- hc$order

    # Filter to non-reference organisms, preserve dendrogram order
    dend_organisms <- ref_organisms[dend_order]
    ordered_organisms <- dend_organisms[dend_organisms != REFERENCE_ORG]

    cat("  Dendrogram order:\n")
    for (org in ordered_organisms) {
        cat("    ", org, "\n", sep = "")
    }

} else {
    stop("Unknown ORDERING_METHOD: ", ORDERING_METHOD)
}

# Final order: reference first, then ordered others
ORGANISM_ORDER <- c(REFERENCE_ORG, ordered_organisms)
cat("\nFinal sequence order:", paste(ORGANISM_ORDER, collapse = "  "), "\n\n")

# ==============================================================================
# PROCESSING
# ==============================================================================

total_positions <- sum(sapply(REGIONS_OF_INTEREST, length))
cat("Positions to visualize:", total_positions, "\n\n")

results <- list()

for (gene in GENE_NAMES) {

    positions <- REGIONS_OF_INTEREST[[gene]]

    if (length(positions) == 0) {
        cat("Skipping", gene, "(no positions)\n")
        next
    }

    cat("[", gene, "] Positions: ", paste(positions, collapse = ", "), "\n", sep = "")

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
    seq_names <- names(aligned)

    # Extract organism codes
    organisms <- sub(pattern = "_.*$", replacement = "", x = seq_names)

    # Reorder sequences according to ORGANISM_ORDER
    new_order <- integer(0)
    for (org in ORGANISM_ORDER) {
        idx <- which(organisms == org)
        if (length(idx) > 0) {
            new_order <- c(new_order, idx[1])
        }
    }

    # Add any organisms not in ORGANISM_ORDER at the end (shouldn't happen)
    missing_idx <- setdiff(seq_along(seq_names), new_order)
    if (length(missing_idx) > 0) {
        cat("  WARNING: Organisms not in order:", organisms[missing_idx], "\n")
        new_order <- c(new_order, missing_idx)
    }

    aligned <- aligned[new_order]
    seq_names <- names(aligned)

    # Find reference sequence in reordered alignment
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
    ref_idx <- ref_idx[1]

    # Build position map
    ref_seq <- strsplit(as.character(aligned[ref_idx]), "")[[1]]
    aln_length <- length(ref_seq)

    ungapped_pos <- 0
    pos_map <- integer(sum(ref_seq != "-"))
    for (i in seq_along(ref_seq)) {
        if (ref_seq[i] != "-") {
            ungapped_pos <- ungapped_pos + 1
            pos_map[ungapped_pos] <- i
        }
    }
    ref_length <- ungapped_pos

    # Process each position
    for (pos in positions) {

        cat("  Position ", pos, " ... ", sep = "")

        if (pos < 1 || pos > ref_length) {
            cat("out of range (max: ", ref_length, ")\n", sep = "")
            results[[length(results) + 1]] <- list(
                gene = gene, position = pos, status = "out_of_range"
            )
            next
        }

        aln_pos <- pos_map[pos]
        residue <- ref_seq[aln_pos]

        start_pos <- max(1, aln_pos - ZOOM_WINDOW)
        end_pos <- min(aln_length, aln_pos + ZOOM_WINDOW)

        cat("aln: ", start_pos, "-", end_pos, " | residue: ", residue, "\n", sep = "")

        # Suppress ggmsa warnings
        p <- suppressWarnings(
            ggmsa::ggmsa(
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
        )

        output_base <- paste0(INPUT_PREFIX, gene, "_pos", pos)

        if (SAVE_PNG) {
            png_path <- file.path(PLOT_OUTPUT_DIR, paste0(output_base, ".png"))
            suppressWarnings(
                ggplot2::ggsave(
                    filename = png_path,
                    plot = p,
                    width = PLOT_WIDTH,
                    height = PLOT_HEIGHT,
                    dpi = 300
                )
            )
        }

        if (SAVE_SVG) {
            svg_path <- file.path(PLOT_OUTPUT_DIR, paste0(output_base, ".svg"))
            suppressWarnings(
                ggplot2::ggsave(
                    filename = svg_path,
                    plot = p,
                    width = PLOT_WIDTH,
                    height = PLOT_HEIGHT
                )
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
