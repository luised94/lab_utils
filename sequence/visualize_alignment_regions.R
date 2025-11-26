# ==============================================================================
# MSA Region Visualization Script
# ==============================================================================
# Purpose: Generate zoomed visualizations of multiple sequence alignments
#          at specified S. cerevisiae positions
#
# Input:   Aligned FASTA files (from alignment generation script)
# Output:  Zoomed region plots (PNG/SVG)
#
# Note:    Positions specified in ungapped S. cerevisiae coordinates
#          Script maps to gapped alignment positions automatically
# Date: 2025-11-26
# ==============================================================================

# === PACKAGE LOADING ==========================================================

library(Biostrings)  # Bioconductor - FASTA I/O, sequence manipulation
library(ggmsa)       # CRAN/Bioconductor - alignment visualization
library(ggplot2)     # CRAN - plot customization

# === FILE PATHS ===============================================================

# Input directory (where aligned FASTA files are)
ALIGNMENT_DIR_path <- "~/data/protein_files/alignments"

# Output directory for plots
PLOT_OUTPUT_DIR_path <- "~/data/protein_files/alignments/plots"

# Create output directory if it doesn't exist
if (!dir.exists(PLOT_OUTPUT_DIR_path)) {
  dir.create(PLOT_OUTPUT_DIR_path, recursive = TRUE)
  cat("Created output directory:", PLOT_OUTPUT_DIR_path, "\n")
}

# === GENES TO PROCESS =========================================================

GENE_NAMES_chr <- c("ORC1", "ORC2", "ORC3", "ORC4", "ORC5", "ORC6",
                    "TOA1", "TOA2")


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

# Window size:  residues around center position
ZOOM_WINDOW_RESIDUES_int <- 10

# === DATASET CONFIGURATION ====================================================

# Which datasets to visualize
VISUALIZE_MODEL_ORGS_lgl <- TRUE
VISUALIZE_UNIREF50_lgl <- TRUE

# UniRef50 subsetting (for cleaner plots)
SUBSET_UNIREF50_lgl <- TRUE
UNIREF50_N_SEQUENCES_int <- 15  # Top N by identity to S. cerevisiae

# === VISUALIZATION SETTINGS ===================================================

# Color scheme for amino acids
# Options: "Chemistry_AA", "Shapely_AA", "Zappo_AA", "Taylor_AA", "Clustal"
PLOT_COLOR_SCHEME_chr <- "Chemistry_AA"

# Show sequence logo (conservation bar) above alignment
SHOW_SEQUENCE_LOGO_lgl <- TRUE

# Show sequence names on left
SHOW_SEQUENCE_NAMES_lgl <- TRUE

# Character width (smaller = more compact)
PLOT_CHAR_WIDTH_dbl <- 0.5

# Plot dimensions
PLOT_WIDTH_inches <- 12
PLOT_HEIGHT_inches <- 8

# Output formats
SAVE_PNG_lgl <- TRUE
SAVE_SVG_lgl <- FALSE

# === S. CEREVISIAE IDENTIFICATION =============================================

# Pattern to identify S. cerevisiae in sequence names
SCER_PATTERN_MODEL_ORGS_chr <- "^Scer_"
SCER_PATTERN_UNIREF50_chr <- "^Sac cer"

# ==============================================================================
# VISUALIZATION WORKFLOW
# ==============================================================================

cat("\n=== STARTING ALIGNMENT VISUALIZATION ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Track visualization outputs
plot_summary_df <- data.frame(
  gene = character(),
  dataset = character(),
  position = integer(),
  status = character(),
  output_file = character(),
  stringsAsFactors = FALSE
)

# Process each gene
for (gene_chr in GENE_NAMES_chr) {

  # Get positions of interest for this gene
  positions_int <- REGIONS_OF_INTEREST_lst[[gene_chr]]

  # Skip if no positions specified
  if (length(positions_int) == 0) {
    cat("Skipping", gene_chr, "(no positions specified)\n\n")
    next
  }

  cat("Processing gene:", gene_chr, "\n")
  cat("Positions to visualize:", paste(positions_int, collapse = ", "), "\n\n")

  # === MODEL ORGANISMS VISUALIZATION ===
  if (VISUALIZE_MODEL_ORGS_lgl) {
    cat("  Dataset: Model Organisms\n")

    # Load alignment
    alignment_path <- file.path(ALIGNMENT_DIR_path,
                                paste0(gene_chr, "_model_orgs_aligned.fasta"))

    if (!file.exists(alignment_path)) {
      cat("    WARNING: Alignment file not found, skipping\n\n")

      for (pos_int in positions_int) {
        plot_summary_df <- rbind(plot_summary_df,
                                 data.frame(gene = gene_chr,
                                          dataset = "model_orgs",
                                          position = pos_int,
                                          status = "alignment_missing",
                                          output_file = NA))
      }

    } else {

      # Load sequences
      aligned_AAStringSet <- Biostrings::readAAStringSet(filepath = alignment_path)

      # Find S. cerevisiae sequence
      seq_names_chr <- names(aligned_AAStringSet)
      scer_index_int <- grep(pattern = SCER_PATTERN_MODEL_ORGS_chr,
                            x = seq_names_chr)

      if (length(scer_index_int) == 0) {
        cat("    WARNING: S. cerevisiae not found in alignment, skipping\n\n")

        for (pos_int in positions_int) {
          plot_summary_df <- rbind(plot_summary_df,
                                   data.frame(gene = gene_chr,
                                            dataset = "model_orgs",
                                            position = pos_int,
                                            status = "scer_not_found",
                                            output_file = NA))
        }

      } else {

        if (length(scer_index_int) > 1) {
          cat("    WARNING: Multiple Scer found, using first\n")
          scer_index_int <- scer_index_int[1]
        }

        # Extract Scer sequence for position mapping
        scer_seq_chr <- as.character(aligned_AAStringSet[scer_index_int])
        scer_seq_split_chr <- strsplit(scer_seq_chr, "")[[1]]

        # Process each position
        for (pos_int in positions_int) {

          cat("    Position", pos_int, "...")

          # Map ungapped Scer position to gapped alignment position
          ungapped_counter_int <- 0
          alignment_position_int <- NA

          for (i in 1:length(scer_seq_split_chr)) {
            if (scer_seq_split_chr[i] != "-") {
              ungapped_counter_int <- ungapped_counter_int + 1
              if (ungapped_counter_int == pos_int) {
                alignment_position_int <- i
                break
              }
            }
          }

          # Check if position was found
          if (is.na(alignment_position_int)) {
            cat(" position out of range, skipping\n")
            plot_summary_df <- rbind(plot_summary_df,
                                     data.frame(gene = gene_chr,
                                              dataset = "model_orgs",
                                              position = pos_int,
                                              status = "position_out_of_range",
                                              output_file = NA))
            next
          }

          # Calculate window boundaries
          start_pos_int <- max(1, alignment_position_int - ZOOM_WINDOW_RESIDUES_int)
          end_pos_int <- min(width(aligned_AAStringSet)[1],
                            alignment_position_int + ZOOM_WINDOW_RESIDUES_int)

          cat(" alignment positions", start_pos_int, "-", end_pos_int, "\n")

          # Generate plot
          plot_obj <- ggmsa::ggmsa(
            msa = aligned_AAStringSet,
            start = start_pos_int,
            end = end_pos_int,
            color = PLOT_COLOR_SCHEME_chr,
            char_width = PLOT_CHAR_WIDTH_dbl,
            seq_name = SHOW_SEQUENCE_NAMES_lgl,
            consensus_views = SHOW_SEQUENCE_LOGO_lgl
          ) +
            ggplot2::ggtitle(paste0(gene_chr, " - Model Organisms - Position ", pos_int)) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"))

          # Save plot
          output_base <- paste0(gene_chr, "_model_orgs_pos", pos_int)

          if (SAVE_PNG_lgl) {
            png_path <- file.path(PLOT_OUTPUT_DIR_path, paste0(output_base, ".png"))
            ggplot2::ggsave(filename = png_path,
                          plot = plot_obj,
                          width = PLOT_WIDTH_inches,
                          height = PLOT_HEIGHT_inches,
                          dpi = 300)
          }

          if (SAVE_SVG_lgl) {
            svg_path <- file.path(PLOT_OUTPUT_DIR_path, paste0(output_base, ".svg"))
            ggplot2::ggsave(filename = svg_path,
                          plot = plot_obj,
                          width = PLOT_WIDTH_inches,
                          height = PLOT_HEIGHT_inches)
          }

          # Record in summary
          plot_summary_df <- rbind(plot_summary_df,
                                   data.frame(gene = gene_chr,
                                            dataset = "model_orgs",
                                            position = pos_int,
                                            status = "complete",
                                            output_file = paste0(output_base, ".png")))
        }
      }
    }

    cat("\n")
  }

  # === UNIREF50 VISUALIZATION ===
  if (VISUALIZE_UNIREF50_lgl) {
    cat("  Dataset: UniRef50\n")

    # Load alignment
    alignment_path <- file.path(ALIGNMENT_DIR_path,
                                paste0(gene_chr, "_uniref50_aligned.fasta"))

    if (!file.exists(alignment_path)) {
      cat("    WARNING: Alignment file not found, skipping\n\n")

      for (pos_int in positions_int) {
        plot_summary_df <- rbind(plot_summary_df,
                                 data.frame(gene = gene_chr,
                                          dataset = "uniref50",
                                          position = pos_int,
                                          status = "alignment_missing",
                                          output_file = NA))
      }

    } else {

      # Load sequences
      aligned_AAStringSet <- Biostrings::readAAStringSet(filepath = alignment_path)

      # Find S. cerevisiae sequence
      seq_names_chr <- names(aligned_AAStringSet)
      scer_index_int <- grep(pattern = SCER_PATTERN_UNIREF50_chr,
                            x = seq_names_chr)

      if (length(scer_index_int) == 0) {
        cat("    WARNING: S. cerevisiae not found in alignment, skipping\n\n")

        for (pos_int in positions_int) {
          plot_summary_df <- rbind(plot_summary_df,
                                   data.frame(gene = gene_chr,
                                            dataset = "uniref50",
                                            position = pos_int,
                                            status = "scer_not_found",
                                            output_file = NA))
        }

      } else {

        if (length(scer_index_int) > 1) {
          cat("    WARNING: Multiple Scer found, using first\n")
          scer_index_int <- scer_index_int[1]
        }

        # Subset sequences if requested
        if (SUBSET_UNIREF50_lgl && length(aligned_AAStringSet) > UNIREF50_N_SEQUENCES_int) {

          cat("    Subsetting to top", UNIREF50_N_SEQUENCES_int, "sequences by identity to Scer\n")

          # Calculate identity to Scer for each sequence
          scer_seq_chr <- as.character(aligned_AAStringSet[scer_index_int])

          identity_scores_dbl <- numeric(length(aligned_AAStringSet))

          for (i in 1:length(aligned_AAStringSet)) {
            seq_chr <- as.character(aligned_AAStringSet[i])

            # Compare position by position
            matches_int <- 0
            valid_positions_int <- 0

            for (j in 1:nchar(seq_chr)) {
              scer_char <- substr(scer_seq_chr, j, j)
              seq_char <- substr(seq_chr, j, j)

              # Count if both are non-gap
              if (scer_char != "-" && seq_char != "-") {
                valid_positions_int <- valid_positions_int + 1
                if (scer_char == seq_char) {
                  matches_int <- matches_int + 1
                }
              }
            }

            # Calculate percent identity
            if (valid_positions_int > 0) {
              identity_scores_dbl[i] <- (matches_int / valid_positions_int) * 100
            } else {
              identity_scores_dbl[i] <- 0
            }
          }

          # Rank by identity (descending)
          ranked_indices_int <- order(identity_scores_dbl, decreasing = TRUE)

          # Keep Scer + top N-1
          if (scer_index_int %in% ranked_indices_int[1:UNIREF50_N_SEQUENCES_int]) {
            # Scer already in top N
            selected_indices_int <- ranked_indices_int[1:UNIREF50_N_SEQUENCES_int]
          } else {
            # Force include Scer
            selected_indices_int <- c(scer_index_int,
                                     ranked_indices_int[1:(UNIREF50_N_SEQUENCES_int - 1)])
          }

          # Subset alignment
          aligned_AAStringSet <- aligned_AAStringSet[selected_indices_int]

          # Update Scer index in subset
          seq_names_chr <- names(aligned_AAStringSet)
          scer_index_int <- grep(pattern = SCER_PATTERN_UNIREF50_chr, x = seq_names_chr)[1]
        }

        # Extract Scer sequence for position mapping
        scer_seq_chr <- as.character(aligned_AAStringSet[scer_index_int])
        scer_seq_split_chr <- strsplit(scer_seq_chr, "")[[1]]

        # Process each position
        for (pos_int in positions_int) {

          cat("    Position", pos_int, "...")

          # Map ungapped Scer position to gapped alignment position
          ungapped_counter_int <- 0
          alignment_position_int <- NA

          for (i in 1:length(scer_seq_split_chr)) {
            if (scer_seq_split_chr[i] != "-") {
              ungapped_counter_int <- ungapped_counter_int + 1
              if (ungapped_counter_int == pos_int) {
                alignment_position_int <- i
                break
              }
            }
          }

          # Check if position was found
          if (is.na(alignment_position_int)) {
            cat(" position out of range, skipping\n")
            plot_summary_df <- rbind(plot_summary_df,
                                     data.frame(gene = gene_chr,
                                              dataset = "uniref50",
                                              position = pos_int,
                                              status = "position_out_of_range",
                                              output_file = NA))
            next
          }

          # Calculate window boundaries
          start_pos_int <- max(1, alignment_position_int - ZOOM_WINDOW_RESIDUES_int)
          end_pos_int <- min(width(aligned_AAStringSet)[1],
                            alignment_position_int + ZOOM_WINDOW_RESIDUES_int)

          cat(" alignment positions", start_pos_int, "-", end_pos_int, "\n")

          # Generate plot
          plot_obj <- ggmsa::ggmsa(
            msa = aligned_AAStringSet,
            start = start_pos_int,
            end = end_pos_int,
            color = PLOT_COLOR_SCHEME_chr,
            char_width = PLOT_CHAR_WIDTH_dbl,
            seq_name = SHOW_SEQUENCE_NAMES_lgl,
            consensus_views = SHOW_SEQUENCE_LOGO_lgl
          ) +
            ggplot2::ggtitle(paste0(gene_chr, " - UniRef50 - Position ", pos_int)) +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"))

          # Save plot
          output_base <- paste0(gene_chr, "_uniref50_pos", pos_int)

          if (SAVE_PNG_lgl) {
            png_path <- file.path(PLOT_OUTPUT_DIR_path, paste0(output_base, ".png"))
            ggplot2::ggsave(filename = png_path,
                          plot = plot_obj,
                          width = PLOT_WIDTH_inches,
                          height = PLOT_HEIGHT_inches,
                          dpi = 300)
          }

          if (SAVE_SVG_lgl) {
            svg_path <- file.path(PLOT_OUTPUT_DIR_path, paste0(output_base, ".svg"))
            ggplot2::ggsave(filename = svg_path,
                          plot = plot_obj,
                          width = PLOT_WIDTH_inches,
                          height = PLOT_HEIGHT_inches)
          }

          # Record in summary
          plot_summary_df <- rbind(plot_summary_df,
                                   data.frame(gene = gene_chr,
                                            dataset = "uniref50",
                                            position = pos_int,
                                            status = "complete",
                                            output_file = paste0(output_base, ".svg")))
        }
      }
    }

    cat("\n")
  }
}

# Save summary
summary_path <- file.path(PLOT_OUTPUT_DIR_path, "visualization_summary.tsv")
write.table(x = plot_summary_df,
            file = summary_path,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

cat("=== VISUALIZATION COMPLETE ===\n")
cat("Summary saved to:", summary_path, "\n")
cat("Total plots generated:", sum(plot_summary_df$status == "complete"), "\n")
print(plot_summary_df)
