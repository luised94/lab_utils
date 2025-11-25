# ==============================================================================
# Script: 03_download_model_organisms.R
# Purpose: Download ORC and TFIIA sequences from model organisms for MSA
# Author: [Your name]
# Date: 2025-01-XX
# ==============================================================================

# CONSTANTS BLOCK ============================================================
# Define all configuration parameters upfront

# Organism taxonomy IDs (17 total)
ORGANISM_TAXIDS_int <- c(
  559292,   # S. cerevisiae (budding yeast)
  27292,    # S. pastorianus (lager yeast)
  1080349,  # S. arboricola
  1080088,  # S. eubayanus
  230603,   # S. uvarum
  27291,    # S. paradoxus
  4954,     # Z. rouxii
  2763761,  # Z. mellis
  37769,    # T. globosa
  5478,     # C. glabrata
  284812,   # S. pombe (fission yeast)
  7227,     # D. melanogaster (fruit fly)
  6239,     # C. elegans (worm)
  10090,    # M. musculus (mouse)
  9606,     # H. sapiens (human)
  8355,     # X. laevis (frog)
  7955      # D. rerio (zebrafish)
)

# Short organism names for FASTA headers
ORGANISM_NAMES_chr <- c(
  "Scer", "Spas", "Sarb", "Seub", "Suva", "Spar",
  "Zrou", "Zmel", "Tglo", "Cgla", "Spom",
  "Dmel", "Cele", "Mmus", "Hsap", "Xlae", "Drer"
)

# Full organism names for metadata
ORGANISM_FULL_NAMES_chr <- c(
  "Saccharomyces cerevisiae",
  "Saccharomyces pastorianus",
  "Saccharomyces arboricola",
  "Saccharomyces eubayanus",
  "Saccharomyces uvarum",
  "Saccharomyces paradoxus",
  "Zygosaccharomyces rouxii",
  "Zygosaccharomyces mellis",
  "Torulaspora globosa",
  "Candida glabrata",
  "Schizosaccharomyces pombe",
  "Drosophila melanogaster",
  "Caenorhabditis elegans",
  "Mus musculus",
  "Homo sapiens",
  "Xenopus laevis",
  "Danio rerio"
)

# Genes to query (8 total)
GENE_NAMES_chr <- c("ORC1", "ORC2", "ORC3", "ORC4", "ORC5", "ORC6", "TOA1", "TOA2")

# UniProt REST API configuration
BASE_URL_uniprot_chr <- "https://rest.uniprot.org"
ENDPOINT_search_chr <- "/uniprotkb/search"
FIELDS_uniprot_chr <- "accession,id,gene_names,organism_name,length,sequence"

# Output configuration
OUTPUT_DIR_path <- "~/data/protein_files"
OUTPUT_PREFIX_chr <- "model_organisms"

# Cache control
FORCE_DOWNLOAD_lgl <- FALSE

# API request configuration
RETRY_MAX_int <- 3
RETRY_DELAY_sec <- 2
REQUEST_DELAY_sec <- 0.5  # Delay between individual queries
USER_AGENT_chr <- "ORC_MSA_Analysis/1.0 (R_script; contact@example.com)"

# LOAD PACKAGES ==============================================================
library(httr2)      # For REST API requests
library(Biostrings) # For FASTA file handling

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
cat("Number of organisms:", length(ORGANISM_TAXIDS_int), "\n")
cat("Number of genes:", length(GENE_NAMES_chr), "\n")
cat("Total queries planned:", length(ORGANISM_TAXIDS_int) * length(GENE_NAMES_chr), "\n")
cat("Output directory:", OUTPUT_DIR_path, "\n")

cat("\nOrganism list:\n")
for (i in 1:length(ORGANISM_NAMES_chr)) {
  cat(sprintf("  %2d. %-5s (%-30s) - %d\n",
              i,
              ORGANISM_NAMES_chr[i],
              ORGANISM_FULL_NAMES_chr[i],
              ORGANISM_TAXIDS_int[i]))
}

cat("\nGene list:\n")
cat(" ", paste(GENE_NAMES_chr, collapse = ", "), "\n")
