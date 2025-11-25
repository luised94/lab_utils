# Load required packages
library(drawProteins)  # Protein domain visualization
library(ggplot2)       # Plot styling and export

# Plot dimensions
FIGURE_WIDTH_inches <- 12
FIGURE_HEIGHT_inches <- 8

# Taxonomy
TAXID_cerevisiae_int <- 559292  # S. cerevisiae

# RELOAD DATA FOR FINAL STYLING ================================================

cat("Reloading curated data for final styling...\n")

CURATED_DOMAINS_path <- "~/data/protein_files/domains_manual_curation_for_plotting.tsv"

domains_curated_df <- read.table(
  CURATED_DOMAINS_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

# Add required columns
gene_order_nls <- c(
  "ORC1" = 1,
  "ORC2" = 2,
  "ORC3" = 3,
  "ORC4" = 4,
  "ORC5" = 5,
  "ORC6" = 6,
  "TOA2" = 7
)

domains_curated_df$order <- gene_order_nls[domains_curated_df$gene_name]

domains_curated_df <- merge(
  domains_curated_df,
  metadata_df[, c("accession", "uniprot_id")],
  by = "accession",
  all.x = TRUE
)

domains_curated_df$entryName <- domains_curated_df$uniprot_id
domains_curated_df$uniprot_id <- NULL
domains_curated_df$taxid <- TAXID_cerevisiae_int
domains_curated_df$source <- "manual_curation"

# Create chains
proteins_in_plot_chr <- unique(domains_curated_df$gene_name)
metadata_for_plot_df <- metadata_df[metadata_df$gene_name %in% proteins_in_plot_chr, ]

chain_entries_df <- data.frame(
  accession = metadata_for_plot_df$accession,
  gene_name = metadata_for_plot_df$gene_name,
  type = "CHAIN",
  description = paste(metadata_for_plot_df$gene_name, "full length"),
  begin = 1,
  end = metadata_for_plot_df$length_aa,
  length = metadata_for_plot_df$length_aa,
  order = gene_order_nls[metadata_for_plot_df$gene_name],
  entryName = metadata_for_plot_df$uniprot_id,
  taxid = TAXID_cerevisiae_int,
  source = "metadata",
  stringsAsFactors = FALSE
)

# Combine
common_cols_chr <- c("accession", "gene_name", "type", "description",
                     "begin", "end", "length", "order", "entryName",
                     "taxid", "source")

plot_data_df <- rbind(
  chain_entries_df[, common_cols_chr],
  domains_curated_df[, common_cols_chr]
)

plot_data_df <- plot_data_df[order(plot_data_df$order, plot_data_df$begin), ]

cat(" Data loaded and prepared\n\n")


# CREATE BASE PLOT =============================================================

cat("Creating base plot layers...\n")

plot_final_plt <- draw_canvas(plot_data_df)

plot_final_plt <- draw_chains(
  p = plot_final_plt,
  data = plot_data_df
)

plot_final_plt <- draw_domains(
  p = plot_final_plt,
  data = plot_data_df,
  label_domains = FALSE
)

plot_final_plt <- draw_motif(
  p = plot_final_plt,
  data = plot_data_df
)

cat(" Base layers created\n\n")


# APPLY FINAL STYLING ==========================================================

cat("Applying final styling...\n")

plot_final_plt <- plot_final_plt +

  # Reverse y-axis so ORC1 appears on top
  scale_y_reverse() +

  # Clean theme
  theme_bw(base_size = 14) +
  theme(
    # Remove all y-axis elements (not informative)
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),

    # Keep x-axis clean
    axis.line.x = element_line(color = "black", size = 0.5),

    # Remove panel elements
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),

    # Legend styling
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key.size = unit(0.8, "cm"),

    # Plot margins
    plot.margin = margin(10, 10, 10, 10)
  ) +

  # Labels
  labs(
    title = "S. cerevisiae ORC Complex - Domain Architecture",
    x = "Amino Acid Position",
    y = NULL,
    fill = "Feature Type"
  )

cat(" Styling applied\n\n")


# ADD PROTEIN NAME LABELS ======================================================

cat("Adding protein name labels at left margin...\n")

# Get y positions for each protein (midpoint of order)
protein_labels_df <- data.frame(
  gene_name = c("ORC1", "ORC2", "ORC3", "ORC4", "ORC5", "ORC6", "TOA2"),
  order = c(1, 2, 3, 4, 5, 6, 7),
  stringsAsFactors = FALSE
)

# Add labels as text annotations
for (i in seq_len(nrow(protein_labels_df))) {
  plot_final_plt <- plot_final_plt +
    annotate(
      "text",
      x = -50,  # Left of plot area
      y = protein_labels_df$order[i],
      label = protein_labels_df$gene_name[i],
      hjust = 1,
      vjust = 0.5,
      size = 5,
      fontface = "bold"
    )
}

cat(" Protein labels added\n\n")


# ADD POSITION LABELS AT CHAIN ENDS ============================================

cat("Adding position labels at chain ends...\n")

# Add start (1) and end (length) labels for each chain
for (i in seq_len(nrow(chain_entries_df))) {
  chain_row <- chain_entries_df[i, ]

  # Start position label
  plot_final_plt <- plot_final_plt +
    annotate(
      "text",
      x = 1,
      y = chain_row$order - 0.35,  # Below chain
      label = "1",
      size = 3,
      color = "gray40"
    )

  # End position label
  plot_final_plt <- plot_final_plt +
    annotate(
      "text",
      x = chain_row$end,
      y = chain_row$order - 0.35,  # Below chain
      label = as.character(chain_row$end),
      size = 3,
      color = "gray40"
    )
}

cat(" Position labels added\n\n")


# SAVE FINAL PUBLICATION-READY PLOT ===========================================

cat("Saving final publication-ready plot...\n")

final_plot_path <- file.path(OUTPUT_DIR_path,
                              paste0(OUTPUT_PREFIX_chr, "_final.pdf"))

ggsave(
  filename = final_plot_path,
  plot = plot_final_plt,
  width = FIGURE_WIDTH_inches,
  height = FIGURE_HEIGHT_inches,
  device = "pdf"
)

cat(" Final plot saved to:", final_plot_path, "\n\n")

# Also save as SVG for easier Illustrator editing
final_svg_path <- file.path(OUTPUT_DIR_path,
                             paste0(OUTPUT_PREFIX_chr, "_final.svg"))

ggsave(
  filename = final_svg_path,
  plot = plot_final_plt,
  width = FIGURE_WIDTH_inches,
  height = FIGURE_HEIGHT_inches,
  device = "svg"
)

cat(" SVG version saved to:", final_svg_path, "\n")
cat("  (Recommended for Illustrator editing)\n\n")


# FINAL VERIFICATION ===========================================================

cat("=" , rep("=", 78), "\n", sep = "")
cat("FINAL PLOT COMPLETE\n")
cat("=" , rep("=", 78), "\n", sep = "")
cat("\nFiles created:\n")
cat("  PDF:", final_plot_path, "\n")
cat("  SVG:", final_svg_path, "\n\n")

cat("Plot features:\n")
cat("   ORC1 on top (reversed y-axis)\n")
cat("   Protein names labeled on left\n")
cat("   Position labels at chain ends (1 and length)\n")
cat("   Y-axis removed (not informative)\n")
cat("   Clean white background\n")
cat("   Domains without internal labels\n")
cat("   Motifs visible\n")
cat("   Legend on right side\n\n")

cat("For Illustrator refinement:\n")
cat("  - Use SVG file for vector editing\n")
cat("  - Add thicker borders to domain rectangles\n")
cat("  - Add emphasis border around motifs (sofr alleles)\n")
cat("  - Adjust colors if needed\n")
cat("  - Fine-tune label positions\n")
cat("  - Add figure legends/annotations\n\n")

cat("Domain and motif counts:\n")
domains_only_df <- plot_data_df[plot_data_df$type == "DOMAIN", ]
motifs_only_df <- plot_data_df[plot_data_df$type == "MOTIF", ]
cat("  Domains:", nrow(domains_only_df), "\n")
cat("  Motifs:", nrow(motifs_only_df), "\n\n")
