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

# BUILD QUERY ================================================================
# Construct query URL with explicit parameters

# Build accession query string (P54784 OR P32833 OR ...)
accession_query_chr <- paste(
  "accession:(",
  paste(ORC_ACCESSIONS_chr, collapse = " OR "),
  ")",
  sep = ""
)

cat("\n=== Building API Query ===\n")
cat("Accession query:", accession_query_chr, "\n")

# Construct full query URL
query_url_chr <- paste0(
  BASE_URL_uniprot_chr,
  ENDPOINT_search_chr,
  "?query=", accession_query_chr,
  "&format=json",
  "&fields=", FIELDS_uniprot_chr
)

cat("Full query URL:\n", query_url_chr, "\n")

# EXECUTE API REQUEST ========================================================
cat("\n=== Executing UniProt API Request ===\n")
cat("Sending GET request...\n")

# Build request with httr2
request_obj <- httr2::request(base_url = query_url_chr) |>
  httr2::req_user_agent(string = USER_AGENT_chr) |>
  httr2::req_retry(
    max_tries = RETRY_MAX_int,
    backoff = ~ RETRY_DELAY_sec
  ) |>
  httr2::req_timeout(seconds = 120)

# Execute request
response_obj <- httr2::req_perform(req = request_obj)

# Check response status
status_code_int <- httr2::resp_status(resp = response_obj)
cat("Response status code:", status_code_int, "\n")

if (status_code_int != 200) {
  stop("API request failed with status code: ", status_code_int)
}

# EXTRACT AND SAVE RAW RESPONSE =============================================
# Get response body as text
response_text_chr <- httr2::resp_body_string(resp = response_obj)

# Save raw JSON response
raw_json_path <- file.path(
  OUTPUT_DIR_path,
  paste0(OUTPUT_PREFIX_chr, "_raw.json")
)

writeLines(text = response_text_chr, con = raw_json_path)

cat("Raw JSON saved to:", raw_json_path, "\n")

# Parse JSON to verify content
response_data_lst <- jsonlite::fromJSON(
  txt = response_text_chr,
  simplifyDataFrame = FALSE
)

# Verify number of results
n_results_int <- length(response_data_lst$results)
cat("Number of proteins retrieved:", n_results_int, "\n")

# Print first protein accession as verification
if (n_results_int > 0) {
  first_accession_chr <- response_data_lst$results[[1]]$primaryAccession
  cat("First protein accession:", first_accession_chr, "\n")
}
