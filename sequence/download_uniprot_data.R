# ==============================================================================
# Script: 01_download_uniprot_data.R
# Purpose: Download S. cerevisiae ORC protein data from UniProt REST API
# Author: [Your name]
# Date: 2025-01-XX
# ==============================================================================

# CONSTANTS BLOCK ============================================================
# Define all configuration parameters upfront
FORCE_DOWNLOAD_lgl <- FALSE

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
  "ft_domain",      # Domain annotations
  "ft_repeat",      # Repeat regions
  "ft_region",      # Regions of interest
  "ft_motif",       # Short sequence motifs
  "ft_chain",       # Mature protein chain
  "xref_interpro",  # InterPro cross-references (for metadata)
  "xref_pfam",      # Pfam cross-references (for metadata)
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

# CHECK CACHE ================================================================
# Skip download if valid cached file exists

raw_json_path <- file.path(
  OUTPUT_DIR_path,
  paste0(OUTPUT_PREFIX_chr, "_raw.json")
)

# Check if cached file exists
if (file.exists(raw_json_path) && !FORCE_DOWNLOAD_lgl) {
  cat("\n=== Found Cached Data ===\n")
  cat("Reading from:", raw_json_path, "\n")

  # Read cached JSON
  response_text_chr <- readLines(con = raw_json_path, warn = FALSE)
  response_text_chr <- paste(response_text_chr, collapse = "\n")

  # Parse to verify
  response_data_lst <- jsonlite::fromJSON(
    txt = response_text_chr,
    simplifyDataFrame = FALSE
  )

  n_results_int <- length(response_data_lst$results)

  # Validate cache
  if (n_results_int == length(ORC_ACCESSIONS_chr)) {
    cat("Cache valid:", n_results_int, "proteins found\n")
    cat("Skipping API request. Set FORCE_DOWNLOAD_lgl = TRUE to refresh.\n")

    # Skip to next chunk (don't execute API request below)
  } else {
    cat("Cache invalid:", n_results_int, "proteins, expected",
        length(ORC_ACCESSIONS_chr), "\n")
    cat("Will re-download...\n")
    FORCE_DOWNLOAD_lgl <- TRUE  # Force download due to invalid cache
  }
} else {
  if (FORCE_DOWNLOAD_lgl) {
    cat("\n=== FORCE_DOWNLOAD enabled - requesting fresh data ===\n")
  } else {
    cat("\n=== No cached data found - downloading ===\n")
  }
}

# Only execute API request if needed
if (!file.exists(raw_json_path) || FORCE_DOWNLOAD_lgl || n_results_int != length(ORC_ACCESSIONS_chr)) {

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

# QUERY INTERPRO API FOR DOMAIN COORDINATES ==================================
cat("\n=== Querying InterPro API for Domain Coordinates ===\n")

# InterPro API configuration
BASE_URL_interpro_chr <- "https://www.ebi.ac.uk/interpro/api"
ENDPOINT_entry_protein_chr <- "/entry/InterPro/protein/UniProt"

# Initialize list to collect domain data
interpro_domains_list <- list()

# Query InterPro for each unique protein accession
unique_accessions_chr <- unique(metadata_df$accession)

for (accession_chr in unique_accessions_chr) {

  # Get gene name for this accession
  gene_chr <- metadata_df$gene_name[metadata_df$accession == accession_chr]

  cat("\nQuerying InterPro for", gene_chr, "(", accession_chr, ")...\n")

  # Build InterPro API URL
  interpro_url_chr <- paste0(
    BASE_URL_interpro_chr,
    ENDPOINT_entry_protein_chr,
    "/",
    accession_chr,
    "/"
  )

  # Execute request with error handling
  tryCatch({
    interpro_request <- httr2::request(base_url = interpro_url_chr) |>
      httr2::req_user_agent(string = USER_AGENT_chr) |>
      httr2::req_retry(max_tries = RETRY_MAX_int) |>
      httr2::req_timeout(seconds = 60)

    interpro_response <- httr2::req_perform(req = interpro_request)

    # Check status
    status_int <- httr2::resp_status(resp = interpro_response)
    if (status_int != 200) {
      cat("  HTTP status:", status_int, "\n")
      next
    }

    interpro_text_chr <- httr2::resp_body_string(resp = interpro_response)

    # Parse with simplifyDataFrame=FALSE to preserve structure
    interpro_data_lst <- jsonlite::fromJSON(
      txt = interpro_text_chr,
      simplifyDataFrame = FALSE
    )

    # Check if results exist
    if (is.null(interpro_data_lst$results) ||
        length(interpro_data_lst$results) == 0) {
      cat("  No InterPro entries found\n")
      next
    }

    cat("  Found", length(interpro_data_lst$results), "InterPro entries\n")

    # Parse each InterPro entry (results is a list)
    for (i in 1:length(interpro_data_lst$results)) {
      entry_lst <- interpro_data_lst$results[[i]]

      # Extract entry metadata
      interpro_accession_chr <- entry_lst$metadata$accession
      interpro_name_chr <- entry_lst$metadata$name
      interpro_type_chr <- entry_lst$metadata$type

      # Extract protein matches (proteins is a list of protein objects)
      if (!is.null(entry_lst$proteins) && length(entry_lst$proteins) > 0) {

        # Iterate through proteins list
        for (j in 1:length(entry_lst$proteins)) {
          protein_lst <- entry_lst$proteins[[j]]

          # Check if this is our target protein (case-insensitive)
          if (tolower(protein_lst$accession) == tolower(accession_chr)) {

            # Extract entry_protein_locations (list of location objects)
            if (!is.null(protein_lst$entry_protein_locations) &&
                length(protein_lst$entry_protein_locations) > 0) {

              # Iterate through locations
              for (k in 1:length(protein_lst$entry_protein_locations)) {
                location_lst <- protein_lst$entry_protein_locations[[k]]

                # Extract fragments (list of fragment objects)
                if (!is.null(location_lst$fragments) &&
                    length(location_lst$fragments) > 0) {

                  # Iterate through fragments
                  for (m in 1:length(location_lst$fragments)) {
                    fragment_lst <- location_lst$fragments[[m]]

                    # Create domain record
                    domain_record <- data.frame(
                      accession = accession_chr,
                      gene_name = gene_chr,
                      source = "InterPro",
                      interpro_id = interpro_accession_chr,
                      type = interpro_type_chr,
                      description = interpro_name_chr,
                      begin = as.integer(fragment_lst$start),
                      end = as.integer(fragment_lst$end),
                      stringsAsFactors = FALSE
                    )

                    interpro_domains_list[[length(interpro_domains_list) + 1]] <-
                      domain_record
                  }
                }
              }
            }
          }
        }
      }
    }

    n_domains_this_protein <- sum(sapply(
      interpro_domains_list,
      function(x) x$accession == accession_chr
    ))
    cat("  Retrieved", n_domains_this_protein, "domain locations\n")

    # Be polite to API
    Sys.sleep(time = 1)

  }, error = function(e) {
    cat("  Error querying", accession_chr, ":", conditionMessage(e), "\n")
  })
}

# Check if we got any domains
if (length(interpro_domains_list) == 0) {
  cat("\nWARNING: No InterPro domains retrieved.\n")
  interpro_domains_df <- data.frame(
    accession = character(0),
    gene_name = character(0),
    source = character(0),
    interpro_id = character(0),
    type = character(0),
    description = character(0),
    begin = integer(0),
    end = integer(0),
    length = integer(0),
    stringsAsFactors = FALSE
  )
} else {
  # Combine all domains into dataframe
  interpro_domains_df <- do.call(rbind, interpro_domains_list)

  # Calculate domain length
  interpro_domains_df$length <- interpro_domains_df$end -
    interpro_domains_df$begin + 1

  # Sort by gene name and start position
  interpro_domains_df <- interpro_domains_df[
    order(interpro_domains_df$gene_name, interpro_domains_df$begin),
  ]
}

# Save InterPro domains
interpro_domains_path <- file.path(
  OUTPUT_DIR_path,
  paste0(OUTPUT_PREFIX_chr, "_interpro_domains.tsv")
)

write.table(
  x = interpro_domains_df,
  file = interpro_domains_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\n=== InterPro Domains Saved ===\n")
cat("File:", interpro_domains_path, "\n")

# VERIFICATION ===============================================================
cat("\n=== InterPro Domains Verification ===\n")
cat("Total domain locations:", nrow(interpro_domains_df), "\n")

if (nrow(interpro_domains_df) > 0) {
  cat("\nDomains per protein:\n")
  print(table(interpro_domains_df$gene_name))

  cat("\nDomain types:\n")
  print(table(interpro_domains_df$type))

  cat("\nSample domains (first 20):\n")
  print(head(interpro_domains_df[, c("gene_name", "interpro_id",
                                      "description", "begin", "end", "length")],
             n = 20))

  cat("\nAAA+ domains found:\n")
  aaa_domains <- interpro_domains_df[
    grepl(pattern = "AAA", x = interpro_domains_df$description, ignore.case = TRUE),
  ]
  if (nrow(aaa_domains) > 0) {
    print(aaa_domains[, c("gene_name", "description", "begin", "end")])
  } else {
    cat("  No AAA+ domains found\n")
  }
}

# QUERY TED API FOR STRUCTURE-BASED DOMAINS =================================
cat("\n=== Querying TED API for Structure-Based Domains ===\n")

# TED API configuration
BASE_URL_ted_chr <- "https://ted.cathdb.info/api/v1"
ENDPOINT_uniprot_summary_chr <- "/uniprot/summary"

# Identify proteins needing TED annotation (missing AAA+/P-loop)
proteins_with_aaa_ploop_chr <- unique(c(
  interpro_domains_df$gene_name[grepl("AAA", interpro_domains_df$description, ignore.case = TRUE)],
  interpro_domains_df$gene_name[grepl("P-loop", interpro_domains_df$description, ignore.case = TRUE)]
))

all_proteins_chr <- metadata_df$gene_name
proteins_needing_ted_chr <- setdiff(all_proteins_chr, proteins_with_aaa_ploop_chr)

cat("Proteins needing TED annotation:", paste(proteins_needing_ted_chr, collapse = ", "), "\n")

# Initialize list for TED domains
ted_domains_list <- list()

# Query TED for each protein needing annotation
for (gene_chr in proteins_needing_ted_chr) {

  # Get accession for this gene
  accession_chr <- metadata_df$accession[metadata_df$gene_name == gene_chr]

  cat("\nQuerying TED for", gene_chr, "(", accession_chr, ")...\n")

  # Build TED API URL
  ted_url_chr <- paste0(
    BASE_URL_ted_chr,
    ENDPOINT_uniprot_summary_chr,
    "/",
    accession_chr
  )

  # Execute request with error handling
  tryCatch({
    ted_request <- httr2::request(base_url = ted_url_chr) |>
      httr2::req_user_agent(string = USER_AGENT_chr) |>
      httr2::req_retry(max_tries = RETRY_MAX_int) |>
      httr2::req_timeout(seconds = 60)

    ted_response <- httr2::req_perform(req = ted_request)

    # Check status
    status_int <- httr2::resp_status(resp = ted_response)
    if (status_int != 200) {
      cat("  HTTP status:", status_int, "\n")
      next
    }

    # Parse JSON response
    ted_response_lst <- httr2::resp_body_json(resp = ted_response)

    # Access the 'data' array (not the top-level response)
    if (is.null(ted_response_lst$data) || length(ted_response_lst$data) == 0) {
      cat("  No TED domains found\n")
      next
    }

    ted_domains_array <- ted_response_lst$data
    cat("  Found", length(ted_domains_array), "TED domains\n")

    # Parse each TED domain in the data array
    for (i in 1:length(ted_domains_array)) {
      domain_lst <- ted_domains_array[[i]]

      # Extract domain information
      ted_id_chr <- domain_lst$ted_id
      chopping_chr <- domain_lst$chopping
      cath_label_chr <- domain_lst$cath_label
      consensus_level_chr <- domain_lst$consensus_level

      # Check if chopping exists
      if (is.null(chopping_chr) || is.na(chopping_chr)) {
        cat("  Warning: No chopping for domain", ted_id_chr, "\n")
        next
      }

      # Parse chopping string (format: "71-191_208-273" for multi-segment, or "2-103" for single)
      # Split by underscore for multiple segments
      segments_chr <- strsplit(x = chopping_chr, split = "_", fixed = TRUE)[[1]]

      # Process each segment as a separate domain entry
      for (segment_chr in segments_chr) {
        # Parse begin-end
        coords_chr <- strsplit(x = segment_chr, split = "-", fixed = TRUE)[[1]]
        begin_int <- as.integer(coords_chr[1])
        end_int <- as.integer(coords_chr[2])

        # Create domain record
        domain_record <- data.frame(
          accession = accession_chr,
          gene_name = gene_chr,
          source = "TED",
          interpro_id = ted_id_chr,
          type = paste0("CATH_", cath_label_chr),
          description = paste0("TED domain (CATH ", cath_label_chr, ", ",
                             consensus_level_chr, " confidence)"),
          begin = begin_int,
          end = end_int,
          stringsAsFactors = FALSE
        )

        ted_domains_list[[length(ted_domains_list) + 1]] <- domain_record
      }
    }

    n_domains_this_protein <- sum(sapply(
      ted_domains_list,
      function(x) x$accession == accession_chr
    ))
    cat("  Retrieved", n_domains_this_protein, "domain segments\n")

    # Be polite to API
    Sys.sleep(time = 1)

  }, error = function(e) {
    cat("  Error querying", accession_chr, ":", conditionMessage(e), "\n")
  })
}

# MERGE INTERPRO AND TED DOMAINS ============================================
cat("\n=== Merging InterPro and TED Domains ===\n")

if (length(ted_domains_list) > 0) {
  # Combine TED domains
  ted_domains_df <- do.call(rbind, ted_domains_list)

  # Calculate length
  ted_domains_df$length <- ted_domains_df$end - ted_domains_df$begin + 1

  cat("TED domains retrieved:", nrow(ted_domains_df), "\n")

  # Combine with InterPro domains
  all_domains_df <- rbind(interpro_domains_df, ted_domains_df)

} else {
  cat("No TED domains retrieved, using InterPro only\n")
  all_domains_df <- interpro_domains_df
}

# Sort by gene name and start position
all_domains_df <- all_domains_df[
  order(all_domains_df$gene_name, all_domains_df$begin),
]

# Save combined domains
all_domains_path <- file.path(
  OUTPUT_DIR_path,
  paste0(OUTPUT_PREFIX_chr, "_all_domains.tsv")
)

write.table(
  x = all_domains_df,
  file = all_domains_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\nCombined domains saved to:", all_domains_path, "\n")

# VERIFICATION ===============================================================
cat("\n=== Combined Domains Verification ===\n")
cat("Total domains:", nrow(all_domains_df), "\n")

cat("\nDomains per protein:\n")
print(table(all_domains_df$gene_name))

cat("\nDomains by source:\n")
print(table(all_domains_df$source))

cat("\nSample combined domains:\n")
print(head(all_domains_df[, c("gene_name", "source", "description", "begin", "end")], n = 15))

cat("\nCoverage check - proteins with domains:\n")
proteins_with_domains_chr <- unique(all_domains_df$gene_name)
print(proteins_with_domains_chr)
missing_proteins_chr <- setdiff(all_proteins_chr, proteins_with_domains_chr)
if (length(missing_proteins_chr) > 0) {
  cat("Missing:", paste(missing_proteins_chr, collapse = ", "), "\n")
} else {
  cat("All proteins have domain annotations!\n")
}
