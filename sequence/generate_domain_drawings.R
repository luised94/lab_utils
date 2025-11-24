# ==============================================================================
# Script 2: Generate Domain Architecture Visualizations
# Purpose: Create publication-quality protein domain drawings for ORC complex
# ==============================================================================

# CONSTANTS ====================================================================

# Input file paths
INPUT_METADATA_path <- "~/data/protein_files/orc_cerevisiae_metadata.tsv"
INPUT_DOMAINS_path <- "~/data/protein_files/orc_cerevisiae_all_domains.tsv"

# Output configuration
OUTPUT_DIR_path <- "~/data/protein_files"
OUTPUT_PREFIX_chr <- "orc_domains"
OUTPUT_FORMATS_chr <- c("svg", "png")

# Plot dimensions
FIGURE_WIDTH_inches <- 12
FIGURE_HEIGHT_inches <- 8

# Domain filtering (will apply after manual consolidation)
MIN_DOMAIN_LENGTH_aa <- 20

# Taxonomy
TAXID_cerevisiae_int <- 559292  # S. cerevisiae


# SETUP ========================================================================

# Load required packages
library(drawProteins)  # Protein domain visualization
library(ggplot2)       # Plot styling and export

cat("Packages loaded successfully\n\n")


# LOAD DATA ====================================================================

# Read metadata file
metadata_df <- read.table(
  INPUT_METADATA_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

cat("Metadata loaded:\n")
print(metadata_df[, c("accession", "gene_name", "uniprot_id", "length_aa")])
cat("\n")

# Read combined domain annotations
domains_df <- read.table(
  INPUT_DOMAINS_path,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

cat("Domain data loaded:\n")
cat("Total domains:", nrow(domains_df), "\n")
cat("Columns:", paste(colnames(domains_df), collapse = ", "), "\n\n")


# VERIFY DATA INTEGRITY ========================================================

# Check that all domain accessions exist in metadata
accessions_in_domains_chr <- unique(domains_df$accession)
accessions_in_metadata_chr <- metadata_df$accession

missing_accessions_lgl <- !accessions_in_domains_chr %in% accessions_in_metadata_chr

if (any(missing_accessions_lgl)) {
  stop("ERROR: Domain data contains accessions not in metadata: ",
       paste(accessions_in_domains_chr[missing_accessions_lgl], collapse = ", "))
} else {
  cat(" All domain accessions match metadata\n")
}

# Display domain counts per protein
domains_per_protein_df <- as.data.frame(table(domains_df$gene_name))
colnames(domains_per_protein_df) <- c("gene_name", "domain_count")
domains_per_protein_df <- merge(
  domains_per_protein_df,
  metadata_df[, c("gene_name", "length_aa")],
  by = "gene_name"
)

cat("\nDomain counts per protein:\n")
print(domains_per_protein_df[order(domains_per_protein_df$gene_name), ])
cat("\n")

# Display unique domain sources
cat("Domain sources:\n")
print(table(domains_df$source))
cat("\n")

# Show full domain list for manual review
cat("Full domain list (for consolidation review):\n")
cat("=" , rep("=", 78), "\n", sep = "")
domain_summary_df <- domains_df[, c("gene_name", "source", "type", 
                                     "description", "begin", "end", "length")]
domain_summary_df <- domain_summary_df[order(domain_summary_df$gene_name, 
                                               domain_summary_df$begin), ]
print(domain_summary_df, row.names = FALSE)
