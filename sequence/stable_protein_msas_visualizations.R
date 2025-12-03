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
# Options: specific gene name (e.g., "ORC1") or "auto" (picks gene with most sequences)
ORDERING_REFERENCE_GENE <- "auto"

# Genes to exclude from auto-selection (e.g., highly conserved genes)
ORDERING_EXCLUDE_GENES <- c("TOA2")

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
# IDENTITY/SIMILARITY OUTPUT SETTINGS
# ==============================================================================

IDENTITY_OUTPUT_FILE <- file.path(PLOT_OUTPUT_DIR, paste0(INPUT_PREFIX, "identity_similarity.tsv"))
SESSION_INFO_FILE <- file.path(PLOT_OUTPUT_DIR, paste0(INPUT_PREFIX, "session_info.txt"))

# ==============================================================================
# SETUP
# ==============================================================================

if (!dir.exists(PLOT_OUTPUT_DIR)) {
    dir.create(path = PLOT_OUTPUT_DIR, recursive = TRUE)
}

# Load BLOSUM62 substitution matrix (must run download_blosum62.R first)
BLOSUM62_PATH <- file.path(path.expand("~/data/protein_files"), "BLOSUM62.rds")
if (!file.exists(BLOSUM62_PATH)) {
    stop("BLOSUM62 matrix not found. Run download_blosum62.R first.")
}
BLOSUM62 <- readRDS(file = BLOSUM62_PATH)

cat("=== STABLE MSA VISUALIZATION ===\n")
cat("Ordering method:", ORDERING_METHOD, "\n\n")

# ==============================================================================
# DETERMINE SEQUENCE ORDER FROM REFERENCE GENE
# ==============================================================================

# Count sequences per gene
cat("Sequence counts per gene:\n")
gene_seq_counts <- integer(length(GENE_NAMES))
names(gene_seq_counts) <- GENE_NAMES

for (i in seq_along(GENE_NAMES)) {
    gene_path <- file.path(ALIGNMENT_DIR, paste0(INPUT_PREFIX, GENE_NAMES[i], "_aligned.fasta"))
    if (file.exists(gene_path)) {
        gene_seq_counts[i] <- length(Biostrings::readAAStringSet(filepath = gene_path))
    } else {
        gene_seq_counts[i] <- 0
    }
    cat("  ", GENE_NAMES[i], ": ", gene_seq_counts[i], "\n", sep = "")
}

# Resolve reference gene
if (ORDERING_REFERENCE_GENE == "auto") {
    # Pick gene with most sequences, excluding specified genes
    eligible_genes <- setdiff(GENE_NAMES, ORDERING_EXCLUDE_GENES)
    eligible_counts <- gene_seq_counts[eligible_genes]
    ORDERING_REFERENCE_GENE <- names(which.max(eligible_counts))
    cat("\nAuto-selected reference gene:", ORDERING_REFERENCE_GENE,
        "(", max(eligible_counts), "sequences )\n")
    cat("Excluded from selection:", paste(ORDERING_EXCLUDE_GENES, collapse = ", "), "\n")
} else {
    cat("\nUsing specified reference gene:", ORDERING_REFERENCE_GENE, "\n")
}

# ==============================================================================
# CALCULATE IDENTITY AND SIMILARITY TO REFERENCE
# ==============================================================================

cat("Calculating identity and similarity to", REFERENCE_ORG, "...\n\n")

# Pre-allocate results list
identity_similarity_results <- list()

for (gene in GENE_NAMES) {

    cat("  [", gene, "] ", sep = "")

    # Load alignment
    alignment_path <- file.path(ALIGNMENT_DIR, paste0(INPUT_PREFIX, gene, "_aligned.fasta"))

    if (!file.exists(alignment_path)) {
        cat("alignment not found: ", alignment_path, "\nskipping\n")
        next
    }

    aligned <- Biostrings::readAAStringSet(filepath = alignment_path)
    seq_names <- names(aligned)
    organisms <- sub(pattern = "_.*$", replacement = "", x = seq_names)

    # Find reference sequence
    ref_idx <- grep(pattern = REFERENCE_PATTERN, x = seq_names)
    if (length(ref_idx) == 0) {
        cat("reference not found, skipping\n")
        next
    }
    ref_idx <- ref_idx[1]

    # Convert to character matrix
    aln_matrix <- as.matrix(aligned)
    ref_seq <- aln_matrix[ref_idx, ]

    # Identify ungapped reference positions
    ungapped_cols <- which(ref_seq != "-")
    n_ref_positions <- length(ungapped_cols)

    # Calculate for each non-reference sequence
    for (i in seq_along(seq_names)) {
        if (i == ref_idx) next

        org <- organisms[i]
        query_seq <- aln_matrix[i, ]

        # Subset to ungapped reference positions
        ref_residues <- ref_seq[ungapped_cols]
        query_residues <- query_seq[ungapped_cols]

        # Positions where query is also non-gap
        both_present <- query_residues != "-"
        n_compared <- sum(both_present)

        if (n_compared == 0) {
            pct_identity <- NA_real_
            pct_similarity <- NA_real_
        } else {
            ref_at_pos <- ref_residues[both_present]
            query_at_pos <- query_residues[both_present]

            # Identity: exact match
            n_identical <- sum(ref_at_pos == query_at_pos)
            pct_identity <- (n_identical / n_compared) * 100

            # Similarity: BLOSUM62 score > 0
            n_similar <- 0
            for (j in seq_len(n_compared)) {
                aa_ref <- ref_at_pos[j]
                aa_query <- query_at_pos[j]

                # Check if both AAs are in BLOSUM62 matrix
                if (aa_ref %in% rownames(BLOSUM62) && aa_query %in% colnames(BLOSUM62)) {
                    if (BLOSUM62[aa_ref, aa_query] > 0) {
                        n_similar <- n_similar + 1
                    }
                }
            }
            pct_similarity <- (n_similar / n_compared) * 100
        }

        # Store result
        identity_similarity_results[[length(identity_similarity_results) + 1]] <- data.frame(
            gene = gene,
            organism = org,
            n_ref_positions = n_ref_positions,
            n_compared = n_compared,
            pct_identity = round(pct_identity, 2),
            pct_similarity = round(pct_similarity, 2),
            stringsAsFactors = FALSE
        )
    }

    cat(length(seq_names) - 1, "sequences processed\n")
}

# Combine into single data frame
identity_similarity_df <- do.call(rbind, identity_similarity_results)

cat("\nIdentity/similarity calculation complete.\n")
cat("Total rows:", nrow(identity_similarity_df), "\n\n")
cat("\n")
stop("breakpoint...")

# --- Load reference gene alignment ---
ref_gene_path <- file.path(
    ALIGNMENT_DIR,
    paste0(INPUT_PREFIX, ORDERING_REFERENCE_GENE, "_aligned.fasta")
)

if (!file.exists(ref_gene_path)) {
    stop("Reference gene alignment not found: ", ref_gene_path)
}

cat("Computing sequence order from", ORDERING_REFERENCE_GENE, "...\n")
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
