# ==============================================================================
# Purpose: Download ORC and TFIIA sequences via single batched UniProt query
# Strategy: One API call for all genes x all organisms, filter locally
# ==============================================================================

# CONFIGURATION BLOCK =========================================================

# Control flags
DRY_RUN_lgl <- TRUE  # Set to FALSE to execute query

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

# JSON PARSING AND INITIAL VERIFICATION =======================================

cat("=== Parsing JSON Response ===\n")

# Extract response body as text
response_text_chr <- httr2::resp_body_string(resp = response_obj)

# Parse JSON with automatic dataframe flattening
response_data <- jsonlite::fromJSON(txt = response_text_chr, simplifyDataFrame = TRUE)

# Extract results array
results_df <- response_data$results

# Verify we got results
if (is.null(results_df) || nrow(results_df) == 0) {
  stop("No results returned from UniProt query", call. = FALSE)
}

cat("Results retrieved:", nrow(results_df), "sequences\n")
cat("Columns in results:", ncol(results_df), "\n")

# Display column names to understand structure
cat("\nColumn names:\n")
print(names(results_df))

# Check for list-columns that need extraction
list_cols_lgl <- sapply(results_df, is.list)
if (any(list_cols_lgl)) {
  cat("\nList-columns requiring extraction:\n")
  print(names(results_df)[list_cols_lgl])
}

cat("\nJSON parsing complete\n\n")

# FIELD EXTRACTION AND DATA ENRICHMENT ========================================

cat("=== Extracting and Flattening Fields ===\n")

# Auto-detect column types among list-columns
dataframe_cols_chr <- names(results_df)[list_cols_lgl & sapply(results_df, is.data.frame)]
non_dataframe_cols_chr <- names(results_df)[list_cols_lgl & sapply(results_df, function(x) !is.data.frame(x))]

cat("Dataframe list-columns:", paste(dataframe_cols_chr, collapse = ", "), "\n")
cat("Non-dataframe list-columns:", paste(non_dataframe_cols_chr, collapse = ", "), "\n\n")

# ASSERTIONS: Validate expected structure before processing
stopifnot(
  "Expected organism to be a dataframe" = "organism" %in% dataframe_cols_chr,
  "Expected sequence to be a dataframe" = "sequence" %in% dataframe_cols_chr,
  "Expected genes to be non-dataframe list" = "genes" %in% non_dataframe_cols_chr
)

# Process dataframe columns - merge with underscore prefix
for (col_name_chr in dataframe_cols_chr) {
  cat("  Flattening:", col_name_chr, "\n")

  # Extract the nested dataframe
  nested_df <- results_df[[col_name_chr]]

  # Prefix column names with parent column name
  prefixed_names_chr <- paste0(col_name_chr, "_", names(nested_df))
  names(nested_df) <- prefixed_names_chr

  # Merge into main dataframe
  results_df <- cbind(results_df, nested_df)

  cat("    Added", ncol(nested_df), "columns\n")
}

# Process non-dataframe list columns
cat("  Processing non-dataframe list-columns...\n")

for (col_name_chr in non_dataframe_cols_chr) {

  if (col_name_chr == "genes") {
    # Custom extraction for genes
    cat("    Custom extraction: genes\n")

    # Extract primary gene name
    results_df$gene_primary <- sapply(results_df$genes, function(gene_list) {
      if (length(gene_list) > 0 && !is.null(gene_list$geneName)) {
        gene_name_value <- gene_list$geneName$value[1]
        if (!is.null(gene_name_value) && length(gene_name_value) > 0) {
          return(gene_name_value)
        }
      }
      return(NA_character_)
    })

    # Extract gene synonyms as comma-separated string
    results_df$gene_synonyms <- sapply(results_df$genes, function(gene_list) {
      if (length(gene_list) > 0 && !is.null(gene_list$synonyms)) {
        synonyms_df <- gene_list$synonyms[[1]]
        if (!is.null(synonyms_df) && nrow(synonyms_df) > 0) {
          synonym_values <- synonyms_df$value
          return(paste(synonym_values, collapse = ", "))
        }
      }
      return(NA_character_)
    })

  } else {
    # Generic fallback: convert to JSON string for unknown complex structures
    cat("    Generic fallback:", col_name_chr, "(converting to JSON string)\n")

    new_col_name <- paste0(col_name_chr, "_json")
    results_df[[new_col_name]] <- sapply(results_df[[col_name_chr]], function(x) {
      jsonlite::toJSON(x = x, auto_unbox = TRUE)
    })
  }
}

# Add computed boolean fields
cat("  Adding computed fields...\n")

# Reviewed status from entryType
results_df$reviewed_lgl <- grepl(pattern = "reviewed", x = results_df$entryType, ignore.case = TRUE)

# Derive short organism names (Genus + species format: Hsap, Scer, etc.)
results_df$organism_short <- sapply(results_df$organism_scientificName, function(sci_name) {
  # Split on whitespace
  parts_chr <- strsplit(x = sci_name, split = "\\s+")[[1]]

  if (length(parts_chr) >= 2) {
    # First letter of genus + first 3 letters of species
    genus_chr <- substr(x = parts_chr[1], start = 1, stop = 1)
    species_chr <- substr(x = parts_chr[2], start = 1, stop = 3)
    short_name <- paste0(toupper(genus_chr), tolower(species_chr))
    return(short_name)
  } else {
    # Fallback for single-word names
    return(substr(x = sci_name, start = 1, stop = 4))
  }
})

cat("    Added reviewed_lgl and organism_short\n")

# Summary of enriched dataframe
cat("\nEnriched dataframe dimensions:", nrow(results_df), "rows x", ncol(results_df), "columns\n")
cat("Column names:\n")
print(names(results_df))

cat("\nField extraction complete\n\n")


# FILTERING AND PRIORITIZATION ================================================

cat("=== Filtering and Prioritizing Sequences ===\n")

# Group by organism and gene, select best sequence per combination
# Priority: reviewed > exact gene match > longest sequence

# Create exact gene match flag for each target gene
results_df$gene_match_exact <- sapply(seq_len(nrow(results_df)), function(i) {
  gene_primary <- results_df$gene_primary[i]
  if (is.na(gene_primary)) return(FALSE)

  # Check if gene_primary matches any of our target genes (case-insensitive)
  any(tolower(gene_primary) == tolower(GENE_NAMES_chr))
})

cat("Exact gene matches:", sum(results_df$gene_match_exact, na.rm = TRUE), "/", nrow(results_df), "\n")

# Split by (organism_taxonId, gene_primary) combination
results_df$group_key <- paste0(results_df$organism_taxonId, "_", results_df$gene_primary)
grouped_list <- split(x = results_df, f = results_df$group_key)

cat("Unique (organism, gene) combinations:", length(grouped_list), "\n")

# Select best sequence from each group
selected_results_list <- lapply(grouped_list, function(group_df) {
  # Sort by priority: reviewed (desc), exact_match (desc), length (desc)
  sorted_indices <- order(
    -group_df$reviewed_lgl,
    -group_df$gene_match_exact,
    -group_df$sequence_length
  )

  # Return top result
  group_df[sorted_indices[1], , drop = FALSE]
})

# Combine selected results back into dataframe
filtered_df <- do.call(rbind, selected_results_list)
rownames(filtered_df) <- NULL

# CRITICAL: Filter to only target genes
cat("Filtering to target genes only...\n")
cat("  Before filtering:", nrow(filtered_df), "sequences\n")

filtered_df <- filtered_df[
  !is.na(filtered_df$gene_primary) &
  tolower(filtered_df$gene_primary) %in% tolower(GENE_NAMES_chr),
]

cat("  After filtering:", nrow(filtered_df), "sequences (target genes only)\n")

# Sort by organism_short then gene_primary for readability
filtered_df <- filtered_df[order(filtered_df$organism_short, tolower(filtered_df$gene_primary)), ]
rownames(filtered_df) <- NULL

cat("\n--- Coverage Analysis ---\n")

# Genes found (case-insensitive unique)
genes_found_unique_chr <- unique(tolower(filtered_df$gene_primary))
genes_found_count <- sum(tolower(GENE_NAMES_chr) %in% genes_found_unique_chr)

cat("Target genes found:", genes_found_count, "/", length(GENE_NAMES_chr), "\n")

# Show which target genes are present
genes_present <- GENE_NAMES_chr[tolower(GENE_NAMES_chr) %in% genes_found_unique_chr]
genes_missing <- GENE_NAMES_chr[!tolower(GENE_NAMES_chr) %in% genes_found_unique_chr]

if (length(genes_present) > 0) {
  cat("  Present:", paste(genes_present, collapse = ", "), "\n")
}
if (length(genes_missing) > 0) {
  cat("  Missing:", paste(genes_missing, collapse = ", "), "\n")
}

# Organisms found
organisms_found_int <- unique(filtered_df$organism_taxonId)
cat("\nTarget organisms found:", length(organisms_found_int), "/", length(ORGANISM_TAXIDS_int), "\n")

# S. cerevisiae verification (taxid 559292)
scer_genes <- unique(tolower(filtered_df$gene_primary[filtered_df$organism_taxonId == 559292]))
cat("\nS. cerevisiae (taxid 559292) coverage:", length(scer_genes), "/ 8 genes\n")
if (length(scer_genes) > 0) {
  scer_found <- GENE_NAMES_chr[tolower(GENE_NAMES_chr) %in% scer_genes]
  scer_missing <- GENE_NAMES_chr[!tolower(GENE_NAMES_chr) %in% scer_genes]
  cat("  Found:", paste(scer_found, collapse = ", "), "\n")
  if (length(scer_missing) > 0) {
    cat("  Missing:", paste(scer_missing, collapse = ", "), "\n")
  }
}

# Missing organisms
missing_taxids_int <- setdiff(ORGANISM_TAXIDS_int, organisms_found_int)
# Organisms with incomplete gene lists
incomplete_orgs <- organisms_found_int[sapply(organisms_found_int, function(taxid) {
  org_genes <- unique(tolower(filtered_df$gene_primary[filtered_df$organism_taxonId == taxid]))
  length(org_genes) < length(GENE_NAMES_chr)
})]

# Reviewed status
cat("\nReviewed status:\n")
print(table(Reviewed = filtered_df$reviewed_lgl))

# Gene x Organism coverage (target genes only, collapsed by case)
cat("\n--- Gene x Organism Coverage (Target Genes) ---\n")

# Normalize gene names to uppercase for display
filtered_df$gene_normalized <- toupper(filtered_df$gene_primary)

coverage_matrix <- table(
  Gene = filtered_df$gene_normalized,
  Organism = filtered_df$organism_short
)
print(coverage_matrix)

# Per-organism gene lists
cat("\n--- Genes Found Per Organism ---\n")
org_gene_summary <- aggregate(
  gene_normalized ~ organism_short + organism_taxonId,
  data = filtered_df,
  FUN = function(x) paste(sort(unique(x)), collapse = ", ")
)
org_gene_summary <- org_gene_summary[order(org_gene_summary$organism_short), ]

for (i in 1:nrow(org_gene_summary)) {
  cat(sprintf("  %-6s (taxid %6d): %s\n",
              org_gene_summary$organism_short[i],
              org_gene_summary$organism_taxonId[i],
              org_gene_summary$gene_normalized[i]))
}

# Export accession table for manual review
cat("\n--- Exporting Accession Table ---\n")

accession_table_df <- filtered_df[, c(
  "gene_normalized",
  "organism_short",
  "organism_scientificName",
  "organism_taxonId",
  "primaryAccession",
  "uniProtkbId",
  "sequence_length",
  "reviewed_lgl"
)]

accession_table_path <- file.path(OUTPUT_DIR_chr, "accession_table.tsv")
write.table(
  x = accession_table_df,
  file = accession_table_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Accession table saved:", accession_table_path, "\n")
cat("  Total sequences:", nrow(accession_table_df), "\n")
cat("  View with: cat", accession_table_path, "\n")

# Generate diagnostic queries for manual verification
cat("\n=== DIAGNOSTIC QUERIES FOR MANUAL VERIFICATION ===\n")

# 1. Query for completely missing organisms
if (length(missing_taxids_int) > 0) {
  cat("\n1) MISSING ORGANISMS (", length(missing_taxids_int), " taxids):\n", sep = "")
  cat("   TaxIDs:", paste(missing_taxids_int, collapse = ", "), "\n")

  missing_org_clause <- paste0("organism_id:", missing_taxids_int, collapse = " OR ")
  backup_query_chr <- paste0(
    "(", gene_clause_combined_chr, " OR ", protein_clause_combined_chr, ") AND (",
    missing_org_clause,
    ") AND (fragment:false)"
  )

  cat("\n   Query to find sequences for missing organisms:\n")
  cat("   ", backup_query_chr, "\n")
} else {
  cat("\n1) All organisms found - no missing organisms\n")
}

# 2. Query for organisms with incomplete gene lists
if (length(incomplete_orgs) > 0) {
  cat("\n2) INCOMPLETE ORGANISMS (", length(incomplete_orgs), " organisms with <8 genes):\n", sep = "")

  for (taxid in incomplete_orgs) {
    org_short <- unique(filtered_df$organism_short[filtered_df$organism_taxonId == taxid])
    org_genes <- unique(tolower(filtered_df$gene_primary[filtered_df$organism_taxonId == taxid]))
    missing_genes <- GENE_NAMES_chr[!tolower(GENE_NAMES_chr) %in% org_genes]

    cat(sprintf("\n   %s (taxid %d): %d/8 genes\n", org_short, taxid, length(org_genes)))
    cat("     Missing:", paste(missing_genes, collapse = ", "), "\n")

    # Generate query for missing genes in this organism
    missing_gene_clauses <- paste0("gene:", missing_genes)
    incomplete_query_chr <- paste0(
      "(", paste(missing_gene_clauses, collapse = " OR "), ") AND ",
      "(organism_id:", taxid, ") AND (fragment:false)"
    )
    cat("     Query:", incomplete_query_chr, "\n")
  }
} else {
  cat("\n2) All found organisms have complete gene sets (8/8)\n")
}

# 3. Accession-based query for verification
cat("\n3) ACCESSION-BASED QUERY (to verify current hits):\n")
accessions_chr <- filtered_df$primaryAccession
accession_query_chr <- paste0(
  "(",
  paste0("accession:", accessions_chr, collapse = " OR "),
  ")"
)

# Save accession query to file
accession_query_path <- file.path(OUTPUT_DIR_chr, "accession_query.txt")
writeLines(text = accession_query_chr, con = accession_query_path)

cat("   Total accessions:", length(accessions_chr), "\n")
cat("   Query saved to:", accession_query_path, "\n")
cat("   First 200 chars:", substr(accession_query_chr, 1, 200), "...\n")

cat("\n=== END DIAGNOSTIC QUERIES ===\n")
cat("\nFiltering complete\n\n")
