# ==============================================================================
# Purpose: Download ORC and TFIIA sequences via single batched UniProt query
# Strategy: One API call for all genes x all organisms, filter locally
# ==============================================================================

# CONFIGURATION BLOCK =========================================================

# Control flags
DRY_RUN_lgl <- FALSE  # Set to FALSE to execute query

# Load required packages
library(httr2)      # For REST API requests
library(Biostrings) # For FASTA file handling

# Target genes (8 total)
GENE_NAMES_chr <- c("ORC1", "ORC2", "ORC3", "ORC4", "ORC5", "ORC6", "TOA1", "TOA2")

# Protein family names for broader matching
# These capture all subunits without specifying numbers
PROTEIN_FAMILIES_chr <- c(
  "origin recognition complex subunit",
  "Transcription initiation factor IIA subunit"
)

# Target organism taxonomy IDs (21 total)
ORGANISM_TAXIDS_int <- c(
  559292,   # S. cerevisiae
  10090,    # M. musculus
  27288,    # N. castellii
  284812,   # S. pombe
  28985,    # K. lactis
  33169,    # A. gossypii
  3702,     # A. thaliana
  4914,     # K. waltii
  4922,     # P. pastoris
  4931,     # S. bayanus
  4934,     # L. kluyveri
  4950,     # T. delbrueckii
  4952,     # Y. lipolytica
  4956,     # Z. rouxii
  4959,     # D. hansenii
  5141,     # N. crassa
  5476,     # C. albicans
  5478,     # C. glabrata
  7227,     # D. melanogaster
  8355,     # X. laevis
  9606      # H. sapiens
)

# UniProt API configuration
BASE_URL_chr <- "https://rest.uniprot.org"
ENDPOINT_chr <- "/uniprotkb/search"
FIELDS_chr <- "accession,id,gene_names,organism_name,organism_id,length,sequence,reviewed"
SIZE_int <- 500  # Max results to retrieve (generous for ~168 potential sequences)

# Request configuration
TIMEOUT_sec <- 120
USER_AGENT_chr <- "ORC_TFIIA_Batch_Download/1.0 (R; contact@example.com)"

# Output configuration
OUTPUT_DIR_chr <- "~/data/protein_files"
METADATA_FILENAME_chr <- "curated_sequences_metadata.tsv"
SUMMARY_FILENAME_chr <- "curated_sequences_summary.tsv"
FASTA_SUFFIX_chr <- "_curated_organisms.fasta"

# Ensure output directory exists
if (!dir.exists(paths = OUTPUT_DIR_chr)) {
  dir.create(path = OUTPUT_DIR_chr, recursive = TRUE)
}

# Validation assertions
stopifnot(
  "GENE_NAMES_chr must not be empty" = length(GENE_NAMES_chr) > 0,
  "ORGANISM_TAXIDS_int must not be empty" = length(ORGANISM_TAXIDS_int) > 0,
  "ORGANISM_TAXIDS_int must be integers" = is.double(ORGANISM_TAXIDS_int),
  "SIZE_int must be positive and less than 500" = SIZE_int > 0 & SIZE_int <= 500,
  "TIMEOUT_sec must be positive" = TIMEOUT_sec > 0,
  "OUTPUT_DIR_chr must be non-empty string" = nchar(OUTPUT_DIR_chr) > 0,
  "Expected results should not exceed SIZE_int" = 
    (length(GENE_NAMES_chr) * length(ORGANISM_TAXIDS_int) * 2) <= SIZE_int
)

cat("=== Configuration Loaded ===\n")
cat("Target genes:", length(GENE_NAMES_chr), "\n")
cat("Target organisms:", length(ORGANISM_TAXIDS_int), "\n")
cat("Expected max results:", length(GENE_NAMES_chr) * length(ORGANISM_TAXIDS_int), "\n")
cat("Dry run mode:", DRY_RUN_lgl, "\n")
cat("Output directory:", OUTPUT_DIR_chr, "\n\n")

# QUERY CONSTRUCTION ==========================================================

cat("=== Constructing Query ===\n")

# Build gene name clause: (gene:ORC1 OR gene:ORC2 OR ...)
gene_clauses_chr <- paste0("gene:", GENE_NAMES_chr)
gene_clause_combined_chr <- paste0("(", paste(gene_clauses_chr, collapse = " OR "), ")")

# Build protein family clause: (protein_name:"..." OR protein_name:"...")
protein_clauses_chr <- paste0("protein_name:\"", PROTEIN_FAMILIES_chr, "\"")
protein_clause_combined_chr <- paste0("(", paste(protein_clauses_chr, collapse = " OR "), ")")

# Combine gene and protein clauses
gene_protein_clause_chr <- paste0("(", gene_clause_combined_chr, " OR ", protein_clause_combined_chr, ")")

# Build organism clause: (organism_id:559292 OR organism_id:10090 OR ...)
organism_clauses_chr <- paste0("organism_id:", ORGANISM_TAXIDS_int)
organism_clause_chr <- paste0("(", paste(organism_clauses_chr, collapse = " OR "), ")")

# Build fragment filter
fragment_clause_chr <- "(fragment:false)"

# Combine all clauses into final query
query_chr <- paste(gene_protein_clause_chr, "AND", organism_clause_chr, "AND", fragment_clause_chr)

cat("Query length:", nchar(query_chr), "characters\n")

if (DRY_RUN_lgl) {
  cat("\n=== DRY RUN: Query String ===\n")
  cat(query_chr, "\n")
  cat("\n=== DRY RUN: Complete - Set DRY_RUN_lgl <- FALSE to execute ===\n")
  stop("Dry run complete", call. = FALSE)
}

cat("Query constructed successfully\n\n")

# API CALL EXECUTION ==========================================================

cat("=== Executing API Request ===\n")

# Define cache file location
RESPONSE_CACHE_FILE_chr <- file.path(OUTPUT_DIR_chr, "uniprot_response_cache.rds")

# Build request object
request_obj <- httr2::request(base_url = paste0(BASE_URL_chr, ENDPOINT_chr)) |>
  httr2::req_url_query(
    query = query_chr,
    format = "json",
    fields = FIELDS_chr,
    size = SIZE_int
  ) |>
  httr2::req_user_agent(string = USER_AGENT_chr) |>
  httr2::req_timeout(seconds = TIMEOUT_sec)

# Check for cached response (session  disk  query)
if (!exists("response_obj")) {
  if (file.exists(RESPONSE_CACHE_FILE_chr)) {
    cat("Loading cached response from disk...\n")
    response_obj <- readRDS(file = RESPONSE_CACHE_FILE_chr)
    cat("Cache loaded:", RESPONSE_CACHE_FILE_chr, "\n")
  } else {
    cat("No cache found - querying UniProt API...\n")
    response_obj <- httr2::req_perform(req = request_obj)

    # Cache to disk for future sessions
    saveRDS(object = response_obj, file = RESPONSE_CACHE_FILE_chr)
    cat("Response cached to disk:", RESPONSE_CACHE_FILE_chr, "\n")
  }
} else {
  cat("Using existing response_obj from current session\n")
}

# Verify response status
response_status_int <- httr2::resp_status(resp = response_obj)
cat("Response status:", response_status_int, "\n")

if (response_status_int != 200) {
  stop("API request failed with status ", response_status_int, call. = FALSE)
}

cat("Request successful\n\n")
