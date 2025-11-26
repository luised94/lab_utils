# ==============================================================================
# Script: download_curated_protein_sequences.R
# Purpose: Download ORC and TFIIA sequences from model organisms for MSA
# Date: 2025-11-26
# ==============================================================================

# CONSTANTS BLOCK ============================================================
# Define all configuration parameters upfront

# Organism taxonomy IDs (17 total)
# Verify taxids at: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi
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
  4931,      # S. bayanus
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
  9606     # H. sapiens
)
stop("breakpoint")
#ORGANISM_TAXIDS_int <- c(
#  559292,   # S. cerevisiae (budding yeast)
#  284812,   # S. pombe (fission yeast)
#  7227,     # D. melanogaster (fruit fly)
#  6239,     # C. elegans (worm)
#  10090,    # M. musculus (mouse)
#  9606,     # H. sapiens (human)
#  8355,     # X. laevis (frog)
#  7955      # D. rerio (zebrafish)
#)

# Short organism names for FASTA headers
# @TODO: These are not used in the query. Can remove.
#ORGANISM_NAMES_chr <- c(
#  "Scer",
#  "Spom",
#  "Dmel",
#  "Cele",
#  "Mmus",
#  "Hsap",
#  "Xlae",
#  "Drer"
#)

# Full organism names for metadata
#ORGANISM_FULL_NAMES_chr <- c(
#  "Saccharomyces cerevisiae",
#  "Schizosaccharomyces pombe",
#  "Drosophila melanogaster",
#  "Caenorhabditis elegans",
#  "Mus musculus",
#  "Homo sapiens",
#  "Xenopus laevis",
#  "Danio rerio"
#)

# Genes to query (8 total)
GENE_NAMES_chr <- c("ORC1", "ORC2", "ORC3", "ORC4", "ORC5", "ORC6", "TOA1", "TOA2")

# UniProt REST API configuration
BASE_URL_uniprot_chr <- "https://rest.uniprot.org"
ENDPOINT_search_chr <- "/uniprotkb/search"
FIELDS_uniprot_chr <- "accession,id,gene_names,organism_name,length,sequence"

# Output configuration
OUTPUT_DIR_path <- "~/data/protein_files"
OUTPUT_PREFIX_chr <- "curated_organisms"

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


# QUERY UNIPROT WITH IMPROVED SEARCH STRATEGIES =============================
cat("\n=== Querying UniProt for Model Organism Sequences ===\n")

# Define protein name mappings for better searching
PROTEIN_NAMES_chr <- c(
  "ORC1" = "origin recognition complex subunit 1",
  "ORC2" = "origin recognition complex subunit 2",
  "ORC3" = "origin recognition complex subunit 3",
  "ORC4" = "origin recognition complex subunit 4",
  "ORC5" = "origin recognition complex subunit 5",
  "ORC6" = "origin recognition complex subunit 6",
  "TOA1" = "transcription initiation factor IIA subunit 1",
  "TOA2" = "transcription initiation factor IIA subunit 2"
)

# Initialize list to store all sequence results
sequences_list <- list()
sequences_count_int <- 0

# Helper function to try a single query
try_query <- function(query_chr, description_chr) {
  tryCatch({
    request_obj <- httr2::request(
      base_url = paste0(BASE_URL_uniprot_chr, ENDPOINT_search_chr)
    ) |>
      httr2::req_url_query(
        query = query_chr,
        format = "json",
        fields = FIELDS_uniprot_chr,
        size = 1
      ) |>
      httr2::req_user_agent(string = USER_AGENT_chr) |>
      httr2::req_retry(max_tries = RETRY_MAX_int, backoff = ~ RETRY_DELAY_sec) |>
      httr2::req_timeout(seconds = 60)

    response_obj <- httr2::req_perform(req = request_obj)

    if (httr2::resp_status(resp = response_obj) == 200) {
      response_text_chr <- httr2::resp_body_string(resp = response_obj)
      response_data_lst <- jsonlite::fromJSON(
        txt = response_text_chr,
        simplifyDataFrame = FALSE
      )

      if (!is.null(response_data_lst$results) &&
          length(response_data_lst$results) > 0) {
        return(list(
          success = TRUE,
          result = response_data_lst$results[[1]],
          method = description_chr
        ))
      }
    }

    return(list(success = FALSE))

  }, error = function(e) {
    return(list(success = FALSE, error = conditionMessage(e)))
  })
}

# Main query loop
for (gene_idx in 1:length(GENE_NAMES_chr)) {
  gene_name_chr <- GENE_NAMES_chr[gene_idx]
  protein_name_chr <- PROTEIN_NAMES_chr[gene_name_chr]

  cat("\n--- Processing", gene_name_chr,
      paste0("(", gene_idx, "/", length(GENE_NAMES_chr), ")"), "---\n")

  for (org_idx in 1:length(ORGANISM_TAXIDS_int)) {
    taxid_int <- ORGANISM_TAXIDS_int[org_idx]
    #org_name_chr <- ORGANISM_NAMES_chr[org_idx]
    #org_full_chr <- ORGANISM_FULL_NAMES_chr[org_idx]

    cat(sprintf("  %s in %s... ", gene_name_chr))#, org_name_chr))

    result_found <- FALSE
    result_data <- NULL
    search_method <- NULL

    # STRATEGY 1: Exact gene match, reviewed only
    if (!result_found) {
      query_chr <- paste0(
        "(gene_exact:", gene_name_chr, ") ",
        "AND (organism_id:", taxid_int, ") ",
        "AND (reviewed:true)"
      )
      attempt <- try_query(query_chr, "gene_exact+reviewed")

      if (attempt$success) {
        result_found <- TRUE
        result_data <- attempt$result
        search_method <- attempt$method
      }
    }

    # STRATEGY 2: Fuzzy gene match, reviewed only
    if (!result_found) {
      query_chr <- paste0(
        "(gene:", gene_name_chr, ") ",
        "AND (organism_id:", taxid_int, ") ",
        "AND (reviewed:true)"
      )
      attempt <- try_query(query_chr, "gene_fuzzy+reviewed")

      if (attempt$success) {
        result_found <- TRUE
        result_data <- attempt$result
        search_method <- attempt$method
      }
    }

    # STRATEGY 3: Protein name search, reviewed only
    if (!result_found) {
      query_chr <- paste0(
        "(protein_name:\"", protein_name_chr, "\") ",
        "AND (organism_id:", taxid_int, ") ",
        "AND (reviewed:true)"
      )
      attempt <- try_query(query_chr, "protein_name+reviewed")

      if (attempt$success) {
        result_found <- TRUE
        result_data <- attempt$result
        search_method <- attempt$method
      }
    }

    # STRATEGY 4: Fuzzy gene match, include unreviewed
    if (!result_found) {
      # Avoid getting small fragments.
      query_chr <- paste0(
        "(gene:", gene_name_chr, ") ",
        "AND (organism_id:", taxid_int, ") ",
        "AND (fragment:false)"
      )
      attempt <- try_query(query_chr, "gene_fuzzy+unreviewed")

      if (attempt$success) {
        result_found <- TRUE
        result_data <- attempt$result
        search_method <- attempt$method
      }
    }

    # STRATEGY 5: Protein name, include unreviewed
    if (!result_found) {
      query_chr <- paste0(
        "(gene:", gene_name_chr, ") ",
        "AND (organism_id:", taxid_int, ") ",
        "AND (fragment:false)"
      )
      attempt <- try_query(query_chr, "protein_name+unreviewed")

      if (attempt$success) {
        result_found <- TRUE
        result_data <- attempt$result
        search_method <- attempt$method
      }
    }

    # Process result or record as missing
    if (result_found) {
      accession_chr <- result_data$primaryAccession
      uniprot_id_chr <- result_data$uniProtkbId
      sequence_chr <- result_data$sequence$value
      length_int <- result_data$sequence$length

      sequences_count_int <- sequences_count_int + 1
      sequences_list[[sequences_count_int]] <- list(
        gene = gene_name_chr,
        #organism_name = org_name_chr,
        #organism_full = org_full_chr,
        taxid = taxid_int,
        accession = accession_chr,
        uniprot_id = uniprot_id_chr,
        sequence = sequence_chr,
        length = length_int,
        status = "found",
        search_method = search_method
      )

      cat("FOUND -", accession_chr,
          paste0("(", length_int, " aa, ", search_method, ")\n"))

    } else {
      cat("NOT FOUND\n")

      sequences_count_int <- sequences_count_int + 1
      sequences_list[[sequences_count_int]] <- list(
        gene = gene_name_chr,
        #organism_name = org_name_chr,
        #organism_full = org_full_chr,
        taxid = taxid_int,
        accession = NA_character_,
        uniprot_id = NA_character_,
        sequence = NA_character_,
        length = NA_integer_,
        status = "missing",
        search_method = NA_character_
      )
    }

    # Rate limiting
    Sys.sleep(time = REQUEST_DELAY_sec)
  }
}

# SUMMARY ====================================================================
cat("\n=== Query Summary ===\n")
cat("Total queries executed:", length(sequences_list), "\n")

status_counts <- table(sapply(sequences_list, function(x) x$status))
cat("Found:", status_counts["found"], "\n")
cat("Missing:", status_counts["missing"], "\n")

cat("\nRetrieval rate:",
    sprintf("%.1f%%", 100 * status_counts["found"] / length(sequences_list)),
    "\n")

# Show search method effectiveness
found_sequences <- sequences_list[sapply(sequences_list, function(x) x$status == "found")]
if (length(found_sequences) > 0) {
  method_counts <- table(sapply(found_sequences, function(x) x$search_method))
  cat("\nSearch methods used:\n")
  for (method in names(method_counts)) {
    cat(sprintf("  %s: %d\n", method, method_counts[method]))
  }
}

# Show which genes/organisms are most problematic
cat("\nMissing by gene:\n")
missing_by_gene <- table(sapply(
  sequences_list[sapply(sequences_list, function(x) x$status == "missing")],
  function(x) x$gene
))
print(missing_by_gene)

# CREATE METADATA TABLE ======================================================
cat("\n=== Creating Metadata Table ===\n")

# Convert list to dataframe
metadata_df <- data.frame(
  gene = character(length(sequences_list)),
  organism_name = character(length(sequences_list)),
  organism_full = character(length(sequences_list)),
  taxid = integer(length(sequences_list)),
  accession = character(length(sequences_list)),
  uniprot_id = character(length(sequences_list)),
  length = integer(length(sequences_list)),
  status = character(length(sequences_list)),
  search_method = character(length(sequences_list)),
  stringsAsFactors = FALSE
)

for (i in 1:length(sequences_list)) {
  seq_lst <- sequences_list[[i]]
  metadata_df$gene[i] <- seq_lst$gene
  metadata_df$organism_name[i] <- seq_lst$organism_name
  metadata_df$organism_full[i] <- seq_lst$organism_full
  metadata_df$taxid[i] <- seq_lst$taxid
  metadata_df$accession[i] <- ifelse(is.na(seq_lst$accession), NA_character_, seq_lst$accession)
  metadata_df$uniprot_id[i] <- ifelse(is.na(seq_lst$uniprot_id), NA_character_, seq_lst$uniprot_id)
  metadata_df$length[i] <- ifelse(is.na(seq_lst$length), NA_integer_, seq_lst$length)
  metadata_df$status[i] <- seq_lst$status
  metadata_df$search_method[i] <- ifelse(is.na(seq_lst$search_method), NA_character_, seq_lst$search_method)
}

# Sort by gene then organism
metadata_df <- metadata_df[order(metadata_df$gene, metadata_df$organism_name), ]

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

# CREATE SUMMARY MATRIX ======================================================
cat("\n=== Creating Summary Matrix ===\n")

# Create presence/absence matrix
summary_matrix <- matrix(
  data = NA_character_,
  nrow = length(GENE_NAMES_chr),
  ncol = length(ORGANISM_NAMES_chr),
  dimnames = list(GENE_NAMES_chr, ORGANISM_NAMES_chr)
)

for (i in 1:nrow(metadata_df)) {
  gene_chr <- metadata_df$gene[i]
  org_chr <- metadata_df$organism_name[i]
  status_chr <- metadata_df$status[i]

  summary_matrix[gene_chr, org_chr] <- ifelse(
    status_chr == "found",
    "yes",
    "no"
  )
}

# Convert to dataframe for saving
summary_df <- as.data.frame(summary_matrix, stringsAsFactors = FALSE)
summary_df <- cbind(Gene = rownames(summary_df), summary_df)
rownames(summary_df) <- NULL

# Save summary matrix
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

cat("Summary matrix saved to:", summary_path, "\n")

# PRINT SUMMARY ==============================================================
cat("\n=== Coverage Summary ===\n")
print(summary_df)

cat("\n=== Organisms by Coverage ===\n")
org_coverage <- colSums(summary_matrix == "yes", na.rm = TRUE)
org_coverage_df <- data.frame(
  Organism = names(org_coverage),
  Found = org_coverage,
  Total = length(GENE_NAMES_chr),
  Percent = sprintf("%.1f%%", 100 * org_coverage / length(GENE_NAMES_chr))
)
org_coverage_df <- org_coverage_df[order(-org_coverage_df$Found), ]
print(org_coverage_df)

# WRITE PER-PROTEIN FASTA FILES ==============================================
cat("\n=== Writing Per-Protein FASTA Files ===\n")

# Filter to only found sequences
found_sequences_list <- sequences_list[sapply(sequences_list, function(x) x$status == "found")]

cat("Total sequences to write:", length(found_sequences_list), "\n\n")

# Process each gene separately
for (gene_name_chr in GENE_NAMES_chr) {

  # Filter sequences for this gene
  gene_sequences <- found_sequences_list[
    sapply(found_sequences_list, function(x) x$gene == gene_name_chr)
  ]

  if (length(gene_sequences) == 0) {
    cat("  ", gene_name_chr, "- No sequences found, skipping\n")
    next
  }

  cat("  ", gene_name_chr, "- Processing", length(gene_sequences), "sequences... ")

  # Extract sequences and create names
  seq_vector_chr <- sapply(gene_sequences, function(x) x$sequence)
  seq_names_chr <- sapply(gene_sequences, function(x) {
    paste0(
      x$organism_name, "_", x$gene,
      "|", x$accession,
      "|", x$taxid,
      "|", x$length, "aa"
    )
  })

  # Create AAStringSet
  sequences_aas <- Biostrings::AAStringSet(x = seq_vector_chr)
  names(sequences_aas) <- seq_names_chr

  # Write FASTA file
  fasta_path <- file.path(
    OUTPUT_DIR_path,
    paste0(gene_name_chr, "_curated_organisms.fasta")
  )

  Biostrings::writeXStringSet(
    x = sequences_aas,
    filepath = fasta_path,
    format = "fasta"
  )

  cat("Written to:", basename(fasta_path), "\n")
}

# VERIFICATION ===============================================================
cat("\n=== FASTA Files Summary ===\n")

for (gene_name_chr in GENE_NAMES_chr) {
  fasta_path <- file.path(
    OUTPUT_DIR_path,
    paste0(gene_name_chr, "_curated_organisms.fasta")
  )

  if (file.exists(fasta_path)) {
    # Read back to verify
    sequences_aas <- Biostrings::readAAStringSet(filepath = fasta_path)
    cat(sprintf("  %-5s: %2d sequences (", 
                gene_name_chr, 
                length(sequences_aas)))
    cat(paste(sapply(strsplit(names(sequences_aas), "_"), `[`, 1), collapse = ", "))
    cat(")\n")
  } else {
    cat(sprintf("  %-5s: File not created (no sequences)\n", gene_name_chr))
  }
}
