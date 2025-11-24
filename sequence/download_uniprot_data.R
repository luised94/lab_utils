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
# Construct query components explicitly

# Build accession query string (P54784 OR P32833 OR ...)
accession_query_chr <- paste(
  "accession:(",
  paste(ORC_ACCESSIONS_chr, collapse = " OR "),
  ")",
  sep = ""
)

cat("\n=== Building API Query ===\n")
cat("Accession query:", accession_query_chr, "\n")

# EXECUTE API REQUEST ========================================================
cat("\n=== Executing UniProt API Request ===\n")
cat("Sending GET request...\n")

# Build request with httr2 - proper method
base_endpoint_chr <- paste0(BASE_URL_uniprot_chr, ENDPOINT_search_chr)

request_obj <- httr2::request(base_url = base_endpoint_chr) |>
  httr2::req_url_query(
    query = accession_query_chr,
    format = "json",
    fields = FIELDS_uniprot_chr,
    .multi = "comma"
  ) |>
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
  first_gene_chr <- response_data_lst$results[[1]]$genes[[1]]$geneName$value
  cat("First protein - Accession:", first_accession_chr, "Gene:", first_gene_chr, "\n")
}


# PARSE METADATA TABLE =======================================================
cat("\n=== Parsing Protein Metadata ===\n")

# Initialize vectors for each column
accession_vct <- character(length = n_results_int)
uniprot_id_vct <- character(length = n_results_int)
gene_name_vct <- character(length = n_results_int)
protein_name_vct <- character(length = n_results_int)
organism_vct <- character(length = n_results_int)
length_aa_vct <- integer(length = n_results_int)
sequence_vct <- character(length = n_results_int)

# Extract data from each protein entry
for (i in 1:n_results_int) {
  protein_lst <- response_data_lst$results[[i]]
  
  # Basic identifiers
  accession_vct[i] <- protein_lst$primaryAccession
  uniprot_id_vct[i] <- protein_lst$uniProtkbId
  
  # Gene name (first gene name from genes array)
  if (length(protein_lst$genes) > 0 && 
      !is.null(protein_lst$genes[[1]]$geneName$value)) {
    gene_name_vct[i] <- protein_lst$genes[[1]]$geneName$value
  } else {
    gene_name_vct[i] <- NA_character_
  }
  
  # Protein description (recommended name)
  if (!is.null(protein_lst$proteinDescription$recommendedName$fullName$value)) {
    protein_name_vct[i] <- protein_lst$proteinDescription$recommendedName$fullName$value
  } else {
    protein_name_vct[i] <- NA_character_
  }
  
  # Organism
  organism_vct[i] <- protein_lst$organism$scientificName
  
  # Sequence information
  length_aa_vct[i] <- protein_lst$sequence$length
  sequence_vct[i] <- protein_lst$sequence$value
}

# Construct metadata dataframe
metadata_df <- data.frame(
  accession = accession_vct,
  uniprot_id = uniprot_id_vct,
  gene_name = gene_name_vct,
  protein_name = protein_name_vct,
  organism = organism_vct,
  length_aa = length_aa_vct,
  sequence = sequence_vct,
  stringsAsFactors = FALSE
)

# Sort by gene name for consistent ordering
metadata_df <- metadata_df[order(metadata_df$gene_name), ]

# Save metadata table
metadata_path <- file.path(
  OUTPUT_DIR_path,
  paste0(OUTPUT_PREFIX_chr, "_metadata.tsv")
)

write.table(
  x = metadata_df,
  file = metadata_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Metadata table saved to:", metadata_path, "\n")

# VERIFICATION ===============================================================
cat("\n=== Metadata Verification ===\n")
cat("Number of rows:", nrow(metadata_df), "\n")
cat("Columns:", paste(names(metadata_df), collapse = ", "), "\n")
cat("\nGene names:\n")
print(metadata_df[, c("gene_name", "accession", "length_aa")])

cat("\nSequence length summary:\n")
print(summary(metadata_df$length_aa))

cat("\nAny NA values in gene_name?", any(is.na(metadata_df$gene_name)), "\n")
cat("Any NA values in sequence?", any(is.na(metadata_df$sequence)), "\n")

# PARSE FEATURE COORDINATES ==================================================
cat("\n=== Parsing Protein Features ===\n")

# Initialize list to collect features from all proteins
features_list <- list()

# Extract features from each protein entry
for (i in 1:n_results_int) {
  protein_lst <- response_data_lst$results[[i]]
  
  # Get protein identifiers for this entry
  accession_chr <- protein_lst$primaryAccession
  gene_chr <- if (length(protein_lst$genes) > 0) {
    protein_lst$genes[[1]]$geneName$value
  } else {
    NA_character_
  }
  
  # Check if features exist
  if (is.null(protein_lst$features) || length(protein_lst$features) == 0) {
    cat("No features found for", accession_chr, "\n")
    next
  }
  
  # Process each feature
  for (j in 1:length(protein_lst$features)) {
    feature_lst <- protein_lst$features[[j]]
    
    # Extract feature information
    feature_type_chr <- feature_lst$type
    feature_desc_chr <- if (!is.null(feature_lst$description)) {
      feature_lst$description
    } else {
      NA_character_
    }
    
    # Extract coordinates
    start_int <- feature_lst$location$start$value
    end_int <- feature_lst$location$end$value
    
    # Create feature record
    feature_record <- data.frame(
      accession = accession_chr,
      gene_name = gene_chr,
      type = feature_type_chr,
      description = feature_desc_chr,
      begin = start_int,
      end = end_int,
      stringsAsFactors = FALSE
    )
    
    features_list[[length(features_list) + 1]] <- feature_record
  }
}

# Combine all features into single dataframe
features_df <- do.call(rbind, features_list)

# Sort by gene name and start position
features_df <- features_df[order(features_df$gene_name, features_df$begin), ]

# Calculate feature length
features_df$length <- features_df$end - features_df$begin + 1

# Save features table
features_path <- file.path(
  OUTPUT_DIR_path,
  paste0(OUTPUT_PREFIX_chr, "_features.tsv")
)

write.table(
  x = features_df,
  file = features_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Features table saved to:", features_path, "\n")

# VERIFICATION ===============================================================
cat("\n=== Features Verification ===\n")
cat("Total features extracted:", nrow(features_df), "\n")
cat("Feature types:\n")
print(table(features_df$type))

cat("\nFeatures per protein:\n")
print(table(features_df$gene_name))

cat("\nSample features:\n")
print(head(features_df[, c("gene_name", "type", "begin", "end", "length")], n = 10))

cat("\nCoordinate validation:\n")
cat("All begin < end?", all(features_df$begin < features_df$end), "\n")
cat("Any negative coordinates?", any(features_df$begin < 0 | features_df$end < 0), "\n")
