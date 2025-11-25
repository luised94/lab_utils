# ==============================================================================
# Script: download_uniprot_refclusters.R
# Purpose: Download UniRef50 cluster members for ORC and TFIIA proteins
# Date: 2025-11-25
# ==============================================================================

# CONSTANTS BLOCK ============================================================

# S. cerevisiae seed accessions (reference proteins)
SEED_ACCESSIONS_chr <- c(
  "P54784",  # ORC1
  "P32833",  # ORC2
  "P54790",  # ORC3
  "P54791",  # ORC4
  "P50874",  # ORC5
  "P38826",  # ORC6
  "P32776",  # TOA1
  "P32774"   # TOA2
)

GENE_NAMES_chr <- c("ORC1", "ORC2", "ORC3", "ORC4", "ORC5", "ORC6", "TOA1", "TOA2")

# UniProtKB API configuration (corrected endpoint)
BASE_URL_uniprotkb_chr <- "https://rest.uniprot.org/uniprotkb/stream"
CLUSTER_IDENTITY_chr <- "50"  # UniRef50

# Output configuration
OUTPUT_DIR_path <- "~/data/protein_files"
OUTPUT_PREFIX_chr <- "uniref50"

# Cache control
FORCE_DOWNLOAD_lgl <- FALSE

# API configuration
RETRY_MAX_int <- 3
RETRY_DELAY_sec <- 2
USER_AGENT_chr <- "ORC_UniRef_Analysis/1.0 (R_script)"

# LOAD PACKAGES ==============================================================
library(httr2)
library(Biostrings)

# VERIFY ENVIRONMENT =========================================================
if (!dir.exists(path = OUTPUT_DIR_path)) {
  dir.create(path = OUTPUT_DIR_path, recursive = TRUE)
}

cat("\n=== Configuration Summary ===\n")
cat("Genes to query:", length(GENE_NAMES_chr), "\n")
cat("Cluster identity:", CLUSTER_IDENTITY_chr, "%\n")
cat("Output directory:", OUTPUT_DIR_path, "\n")
