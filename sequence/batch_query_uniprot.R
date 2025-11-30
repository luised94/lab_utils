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
  "origin recognition complex",
  "transcription factor iia"
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

# BLAST workflow configuration
BLAST_QUERIES_FASTA_chr <- "blast_queries.fasta"
BLAST_QUERY_ACCESSIONS_chr <- "blast_query_accessions.txt"
BLAST_TARGET_ORGANISMS_chr <- "blast_target_organisms.tsv"
UNIPROT_QUERY_FILE_chr <- "uniprot_query_full.txt"
COVERAGE_SUMMARY_FILE_chr <- "uniprot_coverage_summary.txt"
SCER_TAXID_int <- 559292  # S. cerevisiae

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

# Save full query to file for documentation
query_file_path <- file.path(OUTPUT_DIR_chr, UNIPROT_QUERY_FILE_chr)
writeLines(text = query_chr, con = query_file_path)
cat("Query saved to:", query_file_path, "\n")

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

# Check response status
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

# Categorize list-columns by structure type
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

# Add computed fields
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

# Normalize gene names to uppercase for grouping
results_df$gene_normalized <- toupper(results_df$gene_primary)

cat("    Added reviewed_lgl, organism_short, gene_normalized\n")

# DIAGNOSTIC FLAGS AND FILTERING ==============================================

cat("=== Adding Diagnostic Flags ===\n")

# Configuration for diagnostic flags
PROTEIN_FAMILY_ORC_pattern_chr <- "origin recognition complex"
PROTEIN_FAMILY_TFIIA_pattern_chr <- "transcription.*initiation.*factor.*IIA"
FALSE_POSITIVE_PATTERNS_chr <- c("leucine-rich", "mitochondrial*transporter")

# Flag 1: Does gene_primary match target genes?
results_df$gene_matches_target <- sapply(results_df$gene_primary, function(gene) {
  if (is.na(gene)) return(FALSE)
  tolower(gene) %in% tolower(GENE_NAMES_chr)
})

cat("  gene_matches_target:", sum(results_df$gene_matches_target, na.rm = TRUE), "/", nrow(results_df), "\n")

# Flag 2: Which protein family does this belong to?
# Need to check protein name - construct from uniProtkbId or use entryType description
# For now, use a heuristic based on gene name and check if we need to extract protein description
results_df$protein_family <- sapply(seq_len(nrow(results_df)), function(i) {
  gene <- results_df$gene_primary[i]
  # Use gene_normalized to check prefixes
  gene_norm <- results_df$gene_normalized[i]

  if (is.na(gene_norm)) {
    return("unknown")
  }

  # Check gene name for ORC or TOA/TFIIA patterns
  is_orc <- grepl(pattern = "^ORC[0-9]", x = gene_norm, ignore.case = TRUE)
  is_tfiia <- grepl(pattern = "^TOA[0-9]|^TFIIA|^GTF2A", x = gene_norm, ignore.case = TRUE)

  if (is_orc && is_tfiia) {
    return("both")
  } else if (is_orc) {
    return("ORC")
  } else if (is_tfiia) {
    return("TFIIA")
  } else {
    return("unknown")
  }
})

cat("  protein_family distribution:\n")
print(table(results_df$protein_family))

# Flag 3: Query match source (more specific)
results_df$query_match_source <- ifelse(
  results_df$gene_matches_target,
  "gene_name",
  "protein_name"
)

# Flag 4: Needs verification (matched via protein name only)
results_df$needs_verification_lgl <- !results_df$gene_matches_target

cat("  query_match_source:\n")
print(table(results_df$query_match_source))
cat("  needs_verification:", sum(results_df$needs_verification_lgl), "sequences\n")
cat("  query_match_source:\n")
print(table(results_df$query_match_source))

# Flag 5: Potential false positives
results_df$potential_issue_flag <- sapply(results_df$uniProtkbId, function(id) {
  if (is.na(id)) return(NA_character_)

  # Check against false positive patterns
  id_lower <- tolower(id)
  issues <- character(0)

  for (pattern in FALSE_POSITIVE_PATTERNS_chr) {
    if (grepl(pattern = pattern, x = id_lower, ignore.case = TRUE)) {
      issues <- c(issues, pattern)
    }
  }

  if (length(issues) > 0) {
    return(paste(issues, collapse = "; "))
  } else {
    return(NA_character_)
  }
})

flagged_count <- sum(!is.na(results_df$potential_issue_flag))
cat("  potential_issue_flag:", flagged_count, "sequences flagged\n")

if (flagged_count > 0) {
  cat("    Flagged sequences:\n")
  flagged_rows <- results_df[!is.na(results_df$potential_issue_flag),
                              c("gene_primary", "organism_short", "uniProtkbId", "potential_issue_flag")]
  print(flagged_rows, row.names = FALSE)
}

cat("\nDiagnostic flags complete\n\n")


# SORTING FOR MANUAL REVIEW ===================================================

cat("=== Sorting Results for Manual Review ===\n")

# Sort by: organism -> gene -> reviewed -> length
results_df <- results_df[order(
  results_df$organism_taxonId,
  results_df$gene_normalized,
  -results_df$reviewed_lgl,
  -results_df$sequence_length
), ]

rownames(results_df) <- NULL

cat("Total sequences:", nrow(results_df), "\n")
cat("Sorting complete\n\n")

# COVERAGE ANALYSIS ===========================================================

cat("=== Coverage Analysis ===\n")

# Expected sequences per organism
EXPECTED_GENES_PER_ORGANISM_int <- length(GENE_NAMES_chr)

# Count sequences by organism
organism_counts <- table(results_df$organism_taxonId)

cat("Organisms represented:", length(organism_counts), "/", length(ORGANISM_TAXIDS_int), "\n\n")

# Categorize organisms by sequence count
organisms_complete <- names(organism_counts)[organism_counts == EXPECTED_GENES_PER_ORGANISM_int]
organisms_incomplete <- names(organism_counts)[organism_counts < EXPECTED_GENES_PER_ORGANISM_int & organism_counts > 0]
organisms_excess <- names(organism_counts)[organism_counts > EXPECTED_GENES_PER_ORGANISM_int]
organisms_missing <- setdiff(as.character(ORGANISM_TAXIDS_int), names(organism_counts))

# Report complete organisms
if (length(organisms_complete) > 0) {
  cat("COMPLETE organisms (", EXPECTED_GENES_PER_ORGANISM_int, "/", EXPECTED_GENES_PER_ORGANISM_int, " genes): ",
      length(organisms_complete), "\n", sep = "")
  for (taxid in organisms_complete) {
    org_short <- unique(results_df$organism_short[results_df$organism_taxonId == as.numeric(taxid)])
    cat(sprintf("  %s (taxid %s)\n", org_short, taxid))
  }
}

# Report incomplete organisms
if (length(organisms_incomplete) > 0) {
  cat("\nINCOMPLETE organisms (<", EXPECTED_GENES_PER_ORGANISM_int, " genes): ",
      length(organisms_incomplete), "\n", sep = "")
  for (taxid in organisms_incomplete) {
    org_short <- unique(results_df$organism_short[results_df$organism_taxonId == as.numeric(taxid)])
    count <- organism_counts[taxid]
    org_genes <- unique(results_df$gene_normalized[results_df$organism_taxonId == as.numeric(taxid)])
    missing_genes <- setdiff(toupper(GENE_NAMES_chr), org_genes)

    cat(sprintf("  %s (taxid %s): %d/%d genes\n",
                org_short, taxid, count, EXPECTED_GENES_PER_ORGANISM_int))
    cat("    Present:", paste(sort(org_genes), collapse = ", "), "\n")
    if (length(missing_genes) > 0) {
      cat("    Missing:", paste(missing_genes, collapse = ", "), "\n")
    }
  }
}

# Report organisms with excess sequences
if (length(organisms_excess) > 0) {
  cat("\nEXCESS sequences (>", EXPECTED_GENES_PER_ORGANISM_int, " genes - possible isoforms/duplicates): ",
      length(organisms_excess), "\n", sep = "")
  for (taxid in organisms_excess) {
    org_short <- unique(results_df$organism_short[results_df$organism_taxonId == as.numeric(taxid)])
    count <- organism_counts[taxid]
    cat(sprintf("  %s (taxid %s): %d sequences\n", org_short, taxid, count))
  }
}

# Report completely missing organisms
if (length(organisms_missing) > 0) {
  cat("\nMISSING organisms (0 sequences): ", length(organisms_missing), "\n", sep = "")
  cat("  TaxIDs:", paste(organisms_missing, collapse = ", "), "\n")
}

# S. cerevisiae specific verification
SCER_TAXID_int <- 559292
cat("\n--- S. cerevisiae Verification (taxid ", SCER_TAXID_int, ") ---\n", sep = "")

scer_sequences <- results_df[results_df$organism_taxonId == SCER_TAXID_int, ]
scer_gene_count <- nrow(scer_sequences)

cat("Sequences found:", scer_gene_count, "/", EXPECTED_GENES_PER_ORGANISM_int, "\n")

if (scer_gene_count > 0) {
  scer_genes <- unique(scer_sequences$gene_normalized)
  missing_scer_genes <- setdiff(toupper(GENE_NAMES_chr), scer_genes)

  cat("Genes found:", paste(sort(scer_genes), collapse = ", "), "\n")

  if (length(missing_scer_genes) > 0) {
    cat("MISSING genes:", paste(missing_scer_genes, collapse = ", "), "\n")
  }

  # Show accessions for cerevisiae
  cat("\nS. cerevisiae accessions:\n")
  scer_table <- scer_sequences[, c("gene_normalized", "primaryAccession", "sequence_length", "reviewed_lgl")]
  print(scer_table, row.names = FALSE)
} else {
  cat("WARNING: No S. cerevisiae sequences found!\n")
}

# Reviewed status breakdown
cat("\n--- Reviewed Status ---\n")
print(table(Reviewed = results_df$reviewed_lgl))

# Gene matches target breakdown
cat("\n--- Gene Name Matches ---\n")
print(table(MatchesTarget = results_df$gene_matches_target))

# Protein family breakdown
cat("\n--- Protein Family Assignment ---\n")
print(table(Family = results_df$protein_family))

cat("\nCoverage analysis complete\n\n")

# EXPORT AND DIAGNOSTIC QUERIES ===============================================

cat("=== Exporting Results and Generating Diagnostic Queries ===\n")

# Export comprehensive results table
export_columns_chr <- c(
  "gene_normalized",
  "gene_primary",
  "organism_short",
  "organism_scientificName",
  "organism_taxonId",
  "primaryAccession",
  "uniProtkbId",
  "sequence_length",
  "reviewed_lgl",
  "gene_matches_target",
  "protein_family",
  "query_match_source",
  "needs_verification_lgl",
  "potential_issue_flag",
  "gene_synonyms"
)

comprehensive_table_df <- results_df[, export_columns_chr]
comprehensive_table_path <- file.path(OUTPUT_DIR_chr, "uniprot_comprehensive_results.tsv")

write.table(
  x = comprehensive_table_df,
  file = comprehensive_table_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

cat("Comprehensive table:", comprehensive_table_path, "\n")
cat("  Sequences:", nrow(comprehensive_table_df), "\n")

# Export accession-only table
accession_table_df <- results_df[, c(
  "organism_taxonId",
  "organism_short",
  "gene_normalized",
  "primaryAccession",
  "reviewed_lgl",
  "sequence_length"
)]

accession_table_path <- file.path(OUTPUT_DIR_chr, "uniprot_accession_table.tsv")
write.table(
  x = accession_table_df,
  file = accession_table_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Accession table:", accession_table_path, "\n\n")

# Generate diagnostic queries
cat("--- Diagnostic Queries ---\n")

# Query 1: Missing organisms
if (length(organisms_missing) > 0) {
  cat("1) Missing organisms query\n")

  missing_org_clause <- paste0("organism_id:", organisms_missing, collapse = " OR ")
  missing_org_query_chr <- paste0(
    "(", gene_clause_combined_chr, " OR ", protein_clause_combined_chr, ") AND (",
    missing_org_clause,
    ") AND (fragment:false)"
  )

  missing_org_query_path <- file.path(OUTPUT_DIR_chr, "uniprot_query_missing_organisms.txt")
  writeLines(text = missing_org_query_chr, con = missing_org_query_path)
  cat("   Saved:", missing_org_query_path, "\n")
}

# Query 2: Incomplete organisms
if (length(organisms_incomplete) > 0) {
  cat("2) Incomplete organisms queries\n")

  incomplete_queries_chr <- character(0)
  for (taxid in organisms_incomplete) {
    org_short <- unique(results_df$organism_short[results_df$organism_taxonId == as.numeric(taxid)])
    org_genes <- unique(results_df$gene_normalized[results_df$organism_taxonId == as.numeric(taxid)])
    missing_genes <- setdiff(toupper(GENE_NAMES_chr), org_genes)

    if (length(missing_genes) > 0) {
      missing_gene_clauses <- paste0("gene:", missing_genes, collapse = " OR ")
      incomplete_query_chr <- paste0(
        "(", missing_gene_clauses, ") AND ",
        "(organism_id:", taxid, ") AND (fragment:false)"
      )

      query_entry <- paste0("# ", org_short, " (taxid ", taxid, ") - missing: ",
                           paste(missing_genes, collapse = ", "), "\n",
                           incomplete_query_chr, "\n")
      incomplete_queries_chr <- c(incomplete_queries_chr, query_entry)
    }
  }

  incomplete_queries_path <- file.path(OUTPUT_DIR_chr, "queries_incomplete_organisms.txt")
  writeLines(text = incomplete_queries_chr, con = incomplete_queries_path)
  cat("   Saved:", incomplete_queries_path, "\n")
}

# Query 3: Accession-based verification
cat("3) Accession-based query\n")

accessions_chr <- results_df$primaryAccession
accession_query_chr <- paste0(
  "(",
  paste0("accession:", accessions_chr, collapse = " OR "),
  ")"
)

accession_query_path <- file.path(OUTPUT_DIR_chr, "uniprot_query_accessions.txt")
writeLines(text = accession_query_chr, con = accession_query_path)
cat("   Saved:", accession_query_path, "\n")
cat("   Accessions:", length(accessions_chr), "\n")
# BLAST WORKFLOW OUTPUTS ======================================================

cat("\n=== Generating BLAST Workflow Files ===\n")

# Extract S. cerevisiae sequences for BLAST queries
scer_sequences <- results_df[results_df$organism_taxonId == SCER_TAXID_int, ]

if (nrow(scer_sequences) > 0) {
  cat("Extracting S. cerevisiae query sequences...\n")

  # Create FASTA for BLAST queries
  blast_query_fasta_chr <- character(0)
  blast_query_accessions_chr <- character(0)

  for (i in 1:nrow(scer_sequences)) {
    gene <- scer_sequences$gene_normalized[i]
    accession <- scer_sequences$primaryAccession[i]
    sequence <- scer_sequences$sequence_value[i]

    # FASTA format
    fasta_entry <- paste0(">", gene, "|", accession, "|Scer\n", sequence)
    blast_query_fasta_chr <- c(blast_query_fasta_chr, fasta_entry)

    # Accession list
    blast_query_accessions_chr <- c(blast_query_accessions_chr,
                                    paste(gene, accession, sep = "\t"))
  }

  # Write BLAST query FASTA
  blast_fasta_path <- file.path(OUTPUT_DIR_chr, BLAST_QUERIES_FASTA_chr)
  writeLines(text = blast_query_fasta_chr, con = blast_fasta_path)
  cat("  BLAST queries FASTA:", blast_fasta_path, "\n")
  cat("    Sequences:", length(blast_query_fasta_chr), "\n")

  # Write accession list
  blast_acc_path <- file.path(OUTPUT_DIR_chr, BLAST_QUERY_ACCESSIONS_chr)
  writeLines(text = c("gene\taccession", blast_query_accessions_chr),
             con = blast_acc_path)
  cat("  BLAST query accessions:", blast_acc_path, "\n")

} else {
  cat("WARNING: No S. cerevisiae sequences found - cannot generate BLAST queries\n")
}

# Identify organisms needing BLAST
cat("\nIdentifying organisms for BLAST searches...\n")

organisms_for_blast <- c(
  as.numeric(organisms_incomplete),
  as.numeric(organisms_missing)
)

if (length(organisms_for_blast) > 0) {

  blast_target_df <- data.frame(
    taxid = integer(0),
    organism_name = character(0),
    organism_short = character(0),
    category = character(0),
    found_count = integer(0),
    missing_genes = character(0),
    stringsAsFactors = FALSE
  )

  for (taxid in organisms_for_blast) {
    org_sequences <- results_df[results_df$organism_taxonId == taxid, ]

    if (nrow(org_sequences) > 0) {
      # INCOMPLETE
      org_short <- unique(org_sequences$organism_short)[1]
      org_name <- unique(org_sequences$organism_scientificName)[1]
      found_genes <- unique(org_sequences$gene_normalized)
      missing_genes <- setdiff(toupper(GENE_NAMES_chr), found_genes)

      blast_target_df <- rbind(blast_target_df, data.frame(
        taxid = taxid,
        organism_name = org_name,
        organism_short = org_short,
        category = "INCOMPLETE",
        found_count = length(found_genes),
        missing_genes = paste(missing_genes, collapse = ","),
        stringsAsFactors = FALSE
      ))
    } else {
      # MISSING - need to look up organism name from original taxid list
      # For now, just use taxid
      blast_target_df <- rbind(blast_target_df, data.frame(
        taxid = taxid,
        organism_name = NA_character_,
        organism_short = NA_character_,
        category = "MISSING",
        found_count = 0,
        missing_genes = paste(toupper(GENE_NAMES_chr), collapse = ","),
        stringsAsFactors = FALSE
      ))
    }
  }

  # Write BLAST target organisms
  blast_targets_path <- file.path(OUTPUT_DIR_chr, BLAST_TARGET_ORGANISMS_chr)
  write.table(
    x = blast_target_df,
    file = blast_targets_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE,
    na = ""
  )

  cat("  BLAST target organisms:", blast_targets_path, "\n")
  cat("    Organisms:", nrow(blast_target_df), "\n")
  cat("    INCOMPLETE:", sum(blast_target_df$category == "INCOMPLETE"), "\n")
  cat("    MISSING:", sum(blast_target_df$category == "MISSING"), "\n")

} else {
  cat("  No organisms need BLAST searches (all complete)\n")
}

# Create coverage summary file
cat("\nCreating coverage summary...\n")

summary_lines <- c(
  "=== UniProt Query Coverage Summary ===",
  "",
  paste("Query date:", Sys.Date()),
  paste("Total sequences retrieved:", nrow(results_df)),
  paste("Target genes:", length(GENE_NAMES_chr)),
  paste("Target organisms:", length(ORGANISM_TAXIDS_int)),
  "",
  "=== Categories ===",
  "",
  paste("COMPLETE (8/8 genes):", length(organisms_complete)),
  paste("INCOMPLETE (<8 genes):", length(organisms_incomplete)),
  paste("MISSING (0 genes):", length(organisms_missing)),
  paste("EXCESS (>8 genes):", length(organisms_excess)),
  "",
  "=== Match Quality ===",
  "",
  paste("Matched via gene name:", sum(results_df$gene_matches_target)),
  paste("Matched via protein name:", sum(!results_df$gene_matches_target)),
  paste("Needs verification:", sum(results_df$needs_verification_lgl)),
  "",
  "=== BLAST Workflow ===",
  "",
  paste("Organisms requiring BLAST:", length(organisms_for_blast)),
  paste("S. cerevisiae queries available:", nrow(scer_sequences)),
  ""
)

summary_path <- file.path(OUTPUT_DIR_chr, COVERAGE_SUMMARY_FILE_chr)
writeLines(text = summary_lines, con = summary_path)
cat("  Coverage summary:", summary_path, "\n")
cat("\n=== BLAST workflow files generated ===\n")

cat("\n=== Script Complete ===\n")
cat("Output files in:", OUTPUT_DIR_chr, "\n\n")

cat("Main results:\n")
cat("  - comprehensive_results.tsv (all sequences with flags)\n")
cat("  - accession_table.tsv (accessions only)\n")
cat("  - coverage_summary.txt (human-readable summary)\n\n")

cat("Query documentation:\n")
cat("  - uniprot_query_full.txt (complete query used)\n")
cat("  - uniprot_response_cache.rds (cached API response)\n\n")

cat("BLAST workflow:\n")
cat("  - blast_queries.fasta (S. cerevisiae sequences)\n")
cat("  - blast_query_accessions.txt (query accessions)\n")
cat("  - blast_target_organisms.tsv (organisms needing BLAST)\n\n")

cat("Diagnostic queries:\n")
cat("  - query_missing_organisms.txt (if applicable)\n")
cat("  - queries_incomplete_organisms.txt (if applicable)\n")
cat("  - query_accessions.txt (verify current hits)\n\n")

cat("Next steps:\n")
cat("  1. Review comprehensive_results.tsv\n")
cat("  2. Submit BLAST queries for incomplete/missing organisms\n")
cat("  3. Manual curation for ambiguous cases\n")
cat("  4. Merge results into final accession list\n")
