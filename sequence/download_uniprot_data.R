# ==============================================================================
# Script: 01_download_uniprot_data.R
# Purpose: Download S. cerevisiae ORC protein data from UniProt REST API
# Author: [Your name]
# Date: 2025-01-XX
# ==============================================================================

# CONSTANTS BLOCK ============================================================
# Define all configuration parameters upfront

# UniProt API configuration
BASE_URL_uniprot_chr <- "https://rest.uniprot.org"
ENDPOINT_search_chr <- "/uniprotkb/stream"

# S. cerevisiae ORC complex protein accessions
ORC_ACCESSIONS_chr <- c(
  "P54784",  # ORC1
  "P32833",  # ORC2
  "P54790",  # ORC3
  "P54791",  # ORC4
  "P50874",  # ORC5
  "P38826"   # ORC6
)

# Organism taxonomy
TAXONOMY_ID_cerevisiae_int <- 559292  # S. cerevisiae S288c

# Fields to retrieve from UniProt
FIELDS_uniprot_chr <- paste(
  "accession",
  "id",
  "gene_names",
  "protein_name",
  "organism_name",
  "length",
  "sequence",
  "ft_domain",
  "ft_chain",
  "xref_interpro",
  "xref_pfam",
  sep = ","
)

# Output configuration
OUTPUT_DIR_path <- "~/data/protein_files"
OUTPUT_PREFIX_chr <- "orc_cerevisiae"

# API request configuration
RETRY_MAX_int <- 3
RETRY_DELAY_sec <- 2
USER_AGENT_chr <- "ORC_Analysis/1.0 (R_script; contact@example.com)"

# LOAD PACKAGES ==============================================================
library(httr2)      # For REST API requests
library(jsonlite)   # For JSON parsing
library(Biostrings) # For sequence manipulation

# VERIFY ENVIRONMENT =========================================================
# Check output directory exists, create if needed
if (!dir.exists(path = OUTPUT_DIR_path)) {
  dir.create(path = OUTPUT_DIR_path, recursive = TRUE)
  cat("Created output directory:", OUTPUT_DIR_path, "\n")
} else {
  cat("Output directory exists:", OUTPUT_DIR_path, "\n")
}

# Print configuration for verification
cat("\n=== Configuration Summary ===\n")
cat("UniProt base URL:", BASE_URL_uniprot_chr, "\n")
cat("Taxonomy ID:", TAXONOMY_ID_cerevisiae_int, "\n")
cat("Number of proteins:", length(ORC_ACCESSIONS_chr), "\n")
cat("Output directory:", OUTPUT_DIR_path, "\n")
cat("Fields to retrieve:", FIELDS_uniprot_chr, "\n")
