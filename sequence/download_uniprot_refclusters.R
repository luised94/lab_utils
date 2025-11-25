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

# QUERY UNIREF50 CLUSTERS ====================================================
cat("\n=== Querying UniRef50 Clusters ===\n")

# Process each gene
for (gene_idx in 1:length(GENE_NAMES_chr)) {
  gene_name_chr <- GENE_NAMES_chr[gene_idx]
  seed_accession_chr <- SEED_ACCESSIONS_chr[gene_idx]

  cat("\n--- Processing", gene_name_chr,
      paste0("(", gene_idx, "/", length(GENE_NAMES_chr), ")"), "---\n")
  cat("Seed accession:", seed_accession_chr, "\n")

  # Output FASTA path
  fasta_path <- file.path(
    OUTPUT_DIR_path,
    paste0(gene_name_chr, "_uniref50.fasta")
  )

  # Check cache
  if (file.exists(fasta_path) && !FORCE_DOWNLOAD_lgl) {
    cat("Cached file found, skipping download\n")
    sequences_aas <- Biostrings::readAAStringSet(filepath = fasta_path)
    cat("Loaded", length(sequences_aas), "sequences from cache\n")
    next
  }

  # Build cluster ID (assumption: cluster named after seed)
  cluster_id_chr <- paste0("UniRef", CLUSTER_IDENTITY_chr, "_", seed_accession_chr)
  cat("Cluster ID:", cluster_id_chr, "\n")

  # Build query string - search UniProtKB for cluster members
  query_string_chr <- paste0(
    "uniref_cluster_", CLUSTER_IDENTITY_chr, ":", cluster_id_chr
  )

  cat("Query:", query_string_chr, "\n")
  cat("Querying cluster members... ")

  # Execute request
  tryCatch({
    request_obj <- httr2::request(base_url = BASE_URL_uniprotkb_chr) |>
      httr2::req_url_query(
        query = query_string_chr,
        format = "fasta"
      ) |>
      httr2::req_user_agent(string = USER_AGENT_chr) |>
      httr2::req_retry(max_tries = RETRY_MAX_int, backoff = ~ RETRY_DELAY_sec) |>
      httr2::req_timeout(seconds = 300)

    response_obj <- httr2::req_perform(req = request_obj)
    status_int <- httr2::resp_status(resp = response_obj)

    if (status_int == 200) {
      # Get FASTA response
      fasta_text_chr <- httr2::resp_body_string(resp = response_obj)

      # Check for empty response
      if (nchar(fasta_text_chr) < 50) {
        cat("WARNING: Empty or minimal response. Cluster ID may be incorrect.\n")
        cat("The seed", seed_accession_chr, "might belong to a different cluster.\n")
        next
      }

      # Write to temporary file
      temp_fasta_path <- tempfile(fileext = ".fasta")
      writeLines(text = fasta_text_chr, con = temp_fasta_path)

      # Read as AAStringSet
      sequences_aas <- Biostrings::readAAStringSet(filepath = temp_fasta_path)

      # Clean up temp file
      unlink(temp_fasta_path)

      if (length(sequences_aas) == 0) {
        cat("WARNING: No sequences retrieved. Cluster may be empty.\n")
        next
      }

      cat("SUCCESS -", length(sequences_aas), "sequences retrieved\n")

      # Write to final location
      Biostrings::writeXStringSet(
        x = sequences_aas,
        filepath = fasta_path,
        format = "fasta"
      )

      cat("Saved to:", basename(fasta_path), "\n")

    } else {
      cat("HTTP ERROR:", status_int, "\n")
    }

    # Rate limiting
    Sys.sleep(time = 1)

  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
  })
}

# SUMMARY STATISTICS =========================================================
cat("\n=== UniRef50 Cluster Summary ===\n")

summary_df <- data.frame(
  Gene = character(0),
  Cluster_ID = character(0),
  N_Sequences = integer(0),
  stringsAsFactors = FALSE
)

for (gene_idx in 1:length(GENE_NAMES_chr)) {
  gene_name_chr <- GENE_NAMES_chr[gene_idx]
  seed_accession_chr <- SEED_ACCESSIONS_chr[gene_idx]
  cluster_id_chr <- paste0("UniRef", CLUSTER_IDENTITY_chr, "_", seed_accession_chr)

  fasta_path <- file.path(
    OUTPUT_DIR_path,
    paste0(gene_name_chr, "_uniref50.fasta")
  )

  if (file.exists(fasta_path)) {
    sequences_aas <- Biostrings::readAAStringSet(filepath = fasta_path)
    n_seq_int <- length(sequences_aas)

    summary_df <- rbind(
      summary_df,
      data.frame(
        Gene = gene_name_chr,
        Cluster_ID = cluster_id_chr,
        N_Sequences = n_seq_int,
        stringsAsFactors = FALSE
      )
    )
  }
}

print(summary_df)

cat("\nTotal sequences across all clusters:", sum(summary_df$N_Sequences), "\n")

# Save summary
summary_path <- file.path(
  OUTPUT_DIR_path,
  paste0(OUTPUT_PREFIX_chr, "_summary.tsv")
)

write.table(
  x = summary_df,
  file = summary_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("Summary saved to:", summary_path, "\n")

cat("\n=== Script 5 Complete ===\n")
cat("UniRef50 FASTA files ready for MSA analysis\n")
