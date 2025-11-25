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

# TRANSFORM DOMAIN DATA ========================================================

# Step 1: Add protein order based on ORC subunit number
cat("Adding protein order...\n")

# Create order mapping (ORC1=1, ORC2=2, etc.)
gene_order_nls <- c(
  "ORC1" = 1,
  "ORC2" = 2,
  "ORC3" = 3,
  "ORC4" = 4,
  "ORC5" = 5,
  "ORC6" = 6
)

domains_df$order <- gene_order_nls[domains_df$gene_name]

cat("Order assigned:\n")
print(table(domains_df$gene_name, domains_df$order))
cat("\n")


# Step 2: Add entryName from metadata (drawProteins requirement)
cat("Adding entryName column...\n")

domains_df <- merge(
  domains_df,
  metadata_df[, c("accession", "uniprot_id")],
  by = "accession",
  all.x = TRUE
)

# Rename uniprot_id to entryName for drawProteins compatibility
domains_df$entryName <- domains_df$uniprot_id
domains_df$uniprot_id <- NULL

cat("Sample entryNames:\n")
print(unique(domains_df[, c("gene_name", "entryName")]))
cat("\n")


# Step 3: Add taxonomy ID
domains_df$taxid <- TAXID_cerevisiae_int


# Step 4: Handle discontinuous domains (collapse to continuous span)
cat("Handling discontinuous domains...\n")
cat("Original domain count:", nrow(domains_df), "\n")

# For each unique domain (same accession + interpro_id + description),
# collapse to single entry spanning min(begin) to max(end)

# Create grouping key
domains_df$domain_key_chr <- paste(
  domains_df$accession,
  domains_df$interpro_id,
  domains_df$description,
  sep = "___"
)

# Aggregate: take min begin, max end for each domain group
# Keep first occurrence of other columns
domains_collapsed_df <- data.frame(
  accession = character(),
  gene_name = character(),
  source = character(),
  interpro_id = character(),
  type = character(),
  description = character(),
  begin = integer(),
  end = integer(),
  order = integer(),
  entryName = character(),
  taxid = integer(),
  stringsAsFactors = FALSE
)

for (key_chr in unique(domains_df$domain_key_chr)) {
  subset_df <- domains_df[domains_df$domain_key_chr == key_chr, ]
  
  collapsed_row <- subset_df[1, ]  # Keep first row as template
  collapsed_row$begin <- min(subset_df$begin)
  collapsed_row$end <- max(subset_df$end)
  
  domains_collapsed_df <- rbind(domains_collapsed_df, collapsed_row[, colnames(domains_collapsed_df)])
}

# Recalculate length after collapsing
domains_collapsed_df$length <- domains_collapsed_df$end - domains_collapsed_df$begin + 1

cat("After collapsing discontinuous domains:", nrow(domains_collapsed_df), "\n")
cat("\n")


# Step 5: Filter by minimum domain length
cat("Filtering domains by minimum length (", MIN_DOMAIN_LENGTH_aa, " aa)...\n", sep = "")

domains_filtered_df <- domains_collapsed_df[
  domains_collapsed_df$length >= MIN_DOMAIN_LENGTH_aa,
]

cat("Domains after length filter:", nrow(domains_filtered_df), "\n")
cat("Removed:", nrow(domains_collapsed_df) - nrow(domains_filtered_df), "short domains\n")
cat("\n")


# Step 6: Create CHAIN entries for full-length proteins
cat("Creating CHAIN entries for full-length proteins...\n")

chain_entries_df <- data.frame(
  accession = metadata_df$accession,
  gene_name = metadata_df$gene_name,
  source = "metadata",
  interpro_id = NA,
  type = "CHAIN",
  description = paste(metadata_df$gene_name, "full length"),
  begin = 1,
  end = metadata_df$length_aa,
  length = metadata_df$length_aa,
  order = gene_order_nls[metadata_df$gene_name],
  entryName = metadata_df$uniprot_id,
  taxid = TAXID_cerevisiae_int,
  stringsAsFactors = FALSE
)

cat("CHAIN entries created:\n")
print(chain_entries_df[, c("gene_name", "description", "begin", "end")])
cat("\n")


# Step 7: Combine CHAIN and domain entries
domain_features_df <- rbind(chain_entries_df, domains_filtered_df)

# Sort by order (protein) then begin position
domain_features_df <- domain_features_df[
  order(domain_features_df$order, domain_features_df$begin),
]

cat("Final domain features data frame:\n")
cat("Total rows:", nrow(domain_features_df), "(", nrow(chain_entries_df), 
    "chains +", nrow(domains_filtered_df), "domains)\n")
cat("\n")


# VERIFY TRANSFORMED DATA ======================================================

cat("Verification:\n")
cat("=" , rep("=", 78), "\n", sep = "")

# Check all required columns present
required_cols_chr <- c("type", "description", "begin", "end", "order", 
                       "entryName", "taxid")
missing_cols_chr <- required_cols_chr[!required_cols_chr %in% colnames(domain_features_df)]

if (length(missing_cols_chr) > 0) {
  stop("ERROR: Missing required columns: ", paste(missing_cols_chr, collapse = ", "))
} else {
  cat(" All required drawProteins columns present\n")
}

# Check order values
if (all(domain_features_df$order %in% 1:6)) {
  cat(" Order values valid (1-6)\n")
} else {
  stop("ERROR: Invalid order values found")
}

# Check all proteins have CHAIN entry
chains_per_protein_df <- table(domain_features_df$gene_name[domain_features_df$type == "CHAIN"])
if (all(chains_per_protein_df == 1) && length(chains_per_protein_df) == 6) {
  cat(" All 6 proteins have exactly 1 CHAIN entry\n")
} else {
  cat("WARNING: CHAIN entry counts per protein:\n")
  print(chains_per_protein_df)
}

# Show final domain counts per protein
cat("\nFinal domain counts (excluding CHAIN):\n")
domains_only_df <- domain_features_df[domain_features_df$type != "CHAIN", ]
final_counts_df <- as.data.frame(table(domains_only_df$gene_name))
colnames(final_counts_df) <- c("gene_name", "domain_count")
print(final_counts_df)
cat("\n")

# Display final domain list for review
cat("Final domains to be plotted:\n")
cat("=" , rep("=", 78), "\n", sep = "")
review_df <- domain_features_df[, c("gene_name", "type", "description", 
                                     "begin", "end", "length", "source")]
print(review_df, row.names = FALSE)

# CREATE DRAWING CANVAS ========================================================

cat("Initializing drawProteins canvas...\n")

# draw_canvas() requires the full feature data frame
# It sets up the coordinate system for all subsequent drawing
domain_canvas_plt <- draw_canvas(domain_features_df)

cat(" Canvas created\n")
cat("  Canvas class:", class(domain_canvas_plt), "\n")
cat("\n")


# DRAW PROTEIN CHAINS ==========================================================

cat("Drawing protein chains...\n")

# draw_chains() adds horizontal bars representing full-length proteins
# It automatically uses CHAIN type entries
# Each protein will be scaled proportionally to its length
domain_canvas_plt <- draw_chains(
  domain_canvas_plt,
  labels = domain_features_df
)

cat(" Chains added to canvas\n")
cat("\n")


# DISPLAY BASIC PLOT ===========================================================

cat("Displaying basic plot with chains only...\n")

# Add minimal styling for verification
domain_canvas_plt <- domain_canvas_plt +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "Amino Acid Position",
    y = "Protein"
  )

# Display the plot
print(domain_canvas_plt)

cat("\n")
cat("=" , rep("=", 78), "\n", sep = "")
cat("VERIFICATION CHECKPOINT\n")
cat("=" , rep("=", 78), "\n", sep = "")
cat("\nCheck the plot above:\n")
cat("1. You should see 6 horizontal bars (one per protein)\n")
cat("2. Bars should be stacked vertically (ORC1 at top, ORC6 at bottom)\n")
cat("3. Bar widths should be proportional to protein lengths:\n")

# Show expected lengths for reference
length_reference_df <- metadata_df[order(metadata_df$gene_name), 
                                    c("gene_name", "length_aa")]
print(length_reference_df, row.names = FALSE)

cat("\n4. X-axis should show amino acid positions (starts at 0 or 1)\n")
cat("5. Y-axis should show protein numbers or names\n")
cat("\n")
