# Load required packages
library(drawProteins)  # Protein domain visualization
library(ggplot2)       # Plot styling and export

# Taxonomy
TAXID_cerevisiae_int <- 559292  # S. cerevisiae

# LOAD MANUALLY CURATED DOMAINS AND MOTIFS =====================================

cat("Loading manually curated domain file with motifs...\n")

CURATED_DOMAINS_path <- "~/data/protein_files/domains_manual_curation_for_plotting.tsv"

domains_curated_df <- read.table(
  CURATED_DOMAINS_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

cat("Curated data loaded:\n")
cat("Total entries:", nrow(domains_curated_df), "\n")
cat("Domains:", sum(domains_curated_df$type == "DOMAIN"), "\n")
cat("Motifs:", sum(domains_curated_df$type == "MOTIF"), "\n")
cat("Proteins:", length(unique(domains_curated_df$gene_name)), "\n\n")


# ADD REQUIRED DRAWPROTEINS COLUMNS ============================================

cat("Adding required columns for drawProteins...\n")

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

# Add entryName from metadata
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

cat(" Required columns added\n\n")


# CREATE CHAIN ENTRIES =========================================================

cat("Creating CHAIN entries...\n")

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

cat(" CHAIN entries created\n\n")


# COMBINE ALL DATA =============================================================

common_cols_chr <- c("accession", "gene_name", "type", "description", 
                     "begin", "end", "length", "order", "entryName", 
                     "taxid", "source")

plot_data_df <- rbind(
  chain_entries_df[, common_cols_chr],
  domains_curated_df[, common_cols_chr]
)

plot_data_df <- plot_data_df[order(plot_data_df$order, plot_data_df$begin), ]

cat("Final plot data:\n")
cat("Chains:", sum(plot_data_df$type == "CHAIN"), "\n")
cat("Domains:", sum(plot_data_df$type == "DOMAIN"), "\n")
cat("Motifs:", sum(plot_data_df$type == "MOTIF"), "\n\n")


# CREATE PLOT WITH ALL FEATURES ================================================

cat("Creating plot with domains and motifs...\n")

# Initialize canvas
plot_final_plt <- draw_canvas(plot_data_df)

# Draw chains
plot_final_plt <- draw_chains(
  p = plot_final_plt,
  data = plot_data_df
)

# Draw domains (will draw DOMAIN types)
plot_final_plt <- draw_domains(
  p = plot_final_plt,
  data = plot_data_df,
  label_domains = FALSE  # Remove labels as requested
)

# Draw motifs explicitly (will draw MOTIF types)
plot_final_plt <- draw_motif(
  p = plot_final_plt,
  data = plot_data_df
)

cat(" Chains, domains, and motifs drawn\n\n")


# APPLY BASIC STYLING ==========================================================

plot_final_plt <- plot_final_plt +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),  # Remove border
    axis.line.x = element_line(color = "black"),
    legend.position = "right"
  ) +
  labs(
    title = "Protein Domain Architecture",
    x = "Amino Acid Position",
    y = ""
  )


# SAVE PLOT ====================================================================

cat("Saving plot with motifs...\n")

plot_with_motifs_path <- file.path(OUTPUT_DIR_path, "orc_domains_with_motifs.pdf")

ggsave(
  filename = plot_with_motifs_path,
  plot = plot_final_plt,
  width = 12,
  height = 8,
  device = "pdf"
)

cat(" Plot saved to:", plot_with_motifs_path, "\n\n")


# VERIFICATION =================================================================

cat("=" , rep("=", 78), "\n", sep = "")
cat("VERIFICATION\n")
cat("=" , rep("=", 78), "\n", sep = "")
cat("\nPlease open:", plot_with_motifs_path, "\n\n")
cat("Expected:\n")
cat("1. All domains visible (no labels on domains)\n")
cat("2. Motifs visible as small rectangles or markers\n")
cat("3. Panel border removed\n")
cat("4. Check motif positions:\n\n")

motifs_df <- plot_data_df[plot_data_df$type == "MOTIF", ]
print(motifs_df[, c("gene_name", "description", "begin")])

cat("\n")
cat("Next refinements needed:\n")
cat("  - Reverse protein order (ORC1 on top)\n")
cat("  - Customize colors\n")
cat("  - Add position labels at chain ends\n")
cat("\n")
