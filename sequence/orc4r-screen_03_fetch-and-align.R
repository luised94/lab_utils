#!/usr/bin/env Rscript
# orc4r-screen_03_fetch_and_align.R
# Fetch protein sequences, align with DECIPHER, and calculate conservation
# Publication: Martinez-Rodriguez L. et al 2026
#
# Pipeline:
#   Phase 1 - Fetch sequences from UniProt/NCBI using accession table
#   Phase 2 - Align per-gene FASTAs and compute per-position conservation
#
# Input:  stable_protein_accessions.tsv (accession table in working directory)
# Output: ~/data/protein_files/stable_<GENE>.fasta (raw sequences)
#         ~/data/protein_files/alignments/stable_<GENE>_aligned.fasta
#         ~/data/protein_files/alignments/stable_<GENE>_conservation.tsv
#         ~/data/protein_files/alignments/stable_summary.tsv
#
# Dependencies: Biostrings, DECIPHER (Bioconductor); httr2 (CRAN)
#
# Usage: Set working directory to script location, then source or Rscript.
# Installation:
#   install.packages("httr2")
#   BiocManager::install(c("Biostrings", "DECIPHER"))
#   OR
#   renv::install(c("httr2", "Biostrings", "DECIPHER"))
#
# Note: Update CONTACT_EMAIL below before running. NCBI requires a valid
#       contact email for programmatic access.
# Downstream: orc4r-screen_04_msa_visualization.R reads aligned FASTAs.
# Tested with:
#   R 4.4.3
#   Biostrings 2.74.1
#   DECIPHER 3.2.0
#   httr2 1.2.2

# ==============================================================================
# DEPENDENCIES
# ==============================================================================
library(Biostrings)
library(DECIPHER)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Set to directory containing this script and the input TSV
# Adjust if running from a different working directory
INPUT_PATH <- "."

UNIPROT_BASE_URL <- "https://rest.uniprot.org/uniprotkb"
NCBI_EFETCH_URL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

GENE_NAMES <- c("ORC1", "ORC2", "ORC3", "ORC4", "ORC5", "ORC6", "TOA2")
INPUT_PREFIX <- "stable_"
REFERENCE_PATTERN <- "^Scer_"
FASTA_OUTPUT_DIR <- "~/data/protein_files"
ALIGNMENT_OUTPUT_DIR <- "~/data/protein_files/alignments"

INPUT_TSV <- file.path(INPUT_PATH, "orc4r-screen_01_protein-accessions.tsv")
OUTPUT_DIR <- FASTA_OUTPUT_DIR
QUERIES_OUTPUT <- file.path(INPUT_PATH, "stable_fetch_queries.tsv")

# Rate limit settings (conservative for unkeyed access)
REQUEST_DELAY_SEC <- 0.6
REQUEST_DELAY_JITTER <- 0.1
REQUEST_TIMEOUT_SEC <- 30
MAX_RETRIES <- 3
RETRY_BACKOFF_SEC <- 2

# Contact info for NCBI (required for polite access)
USER_AGENT <- "ProteinFetchScript/1.0"
CONTACT_EMAIL <- "luised94@mit.edu"

if (CONTACT_EMAIL == "luised94@mit.edu") {
    message("NOTE: Update CONTACT_EMAIL to your own address before distributing this script.")
}

# ==============================================================================
# REQUIRED INPUTS
# ==============================================================================

if (!file.exists(INPUT_TSV)) {
    stop("Input file not found: ", INPUT_TSV)
}

accessions <- utils::read.delim(
    file = INPUT_TSV,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
)

# --- Validate input columns ---
required_cols <- c("organism", "gene", "accession", "database")
missing_cols <- setdiff(required_cols, colnames(accessions))
if (length(missing_cols) > 0) {
    stop("Accession TSV missing required columns: ", paste(missing_cols, collapse = ", "),
         "\nExpected: ", paste(required_cols, collapse = ", "))
}

cat("Loaded", nrow(accessions), "accessions from", INPUT_TSV, "\n")
cat("Organisms:", length(unique(accessions$organism)), "\n")
cat("Genes:", length(unique(accessions$gene)), "\n")
cat("UniProt:", sum(accessions$database == "uniprot"), "\n")
cat("NCBI:", sum(accessions$database == "ncbi"), "\n\n")

# ==============================================================================
# PREPROCESSING - BUILD URLS
# ==============================================================================

urls <- character(nrow(accessions))

for (i in seq_len(nrow(accessions))) {
    db <- accessions$database[i]
    acc <- accessions$accession[i]

    if (db == "uniprot") {
        urls[i] <- paste0(UNIPROT_BASE_URL, "/", acc, ".fasta")
    } else if (db == "ncbi") {
        urls[i] <- paste0(
            NCBI_EFETCH_URL,
            "?db=protein",
            "&id=", acc,
            "&rettype=fasta",
            "&retmode=text"
        )
    } else {
        urls[i] <- NA_character_
        warning("Unknown database '", db, "' for accession: ", acc)
    }
}

accessions$url <- urls

# ==============================================================================
# OUTPUT - QUERIES FOR MANUAL TESTING
# ==============================================================================

utils::write.table(
    x = accessions,
    file = QUERIES_OUTPUT,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

cat("Wrote query URLs to:", QUERIES_OUTPUT, "\n\n")

# Print sample URLs for quick manual testing
cat("=== SAMPLE URLS FOR MANUAL TESTING ===\n\n")

cat("-- UniProt example --\n")
uniprot_idx <- which(accessions$database == "uniprot")[1]
cat(accessions$organism[uniprot_idx], accessions$gene[uniprot_idx], "\n")
cat(accessions$url[uniprot_idx], "\n\n")

cat("-- NCBI example --\n")
ncbi_idx <- which(accessions$database == "ncbi")[1]
cat(accessions$organism[ncbi_idx], accessions$gene[ncbi_idx], "\n")
cat(accessions$url[ncbi_idx], "\n\n")

# ==============================================================================
# EXECUTION - FETCH SEQUENCES
# ==============================================================================

cat("=== FETCHING SEQUENCES ===\n")
cat("Delay:", REQUEST_DELAY_SEC, "s (+/-", REQUEST_DELAY_JITTER, "s jitter)\n")
cat("Estimated time:", round(nrow(accessions) * REQUEST_DELAY_SEC / 60, 1), "min\n\n")

# Storage for results
sequences <- vector(mode = "list", length = nrow(accessions))
failures <- character(0)

for (i in seq_len(nrow(accessions))) {
    org <- accessions$organism[i]
    gene <- accessions$gene[i]
    acc <- accessions$accession[i]
    db <- accessions$database[i]
    url <- accessions$url[i]

    cat("[", i, "/", nrow(accessions), "] ", org, " ", gene, " (", acc, ")... ", sep = "")

    # Build request with httr2
    req <- httr2::request(base_url = url)
    req <- httr2::req_user_agent(req = req, string = USER_AGENT)
    req <- httr2::req_timeout(req = req, seconds = REQUEST_TIMEOUT_SEC)
    req <- httr2::req_retry(
        req = req,
        max_tries = MAX_RETRIES,
        backoff = ~ RETRY_BACKOFF_SEC
    )

    # Add NCBI-specific headers
    if (db == "ncbi") {
        req <- httr2::req_headers(req, `tool` = USER_AGENT, `email` = CONTACT_EMAIL)
    }

    # Perform request
    resp <- tryCatch(
        expr = httr2::req_perform(req = req),
        error = function(e) {
            cat("ERROR:", conditionMessage(e), "\n")
            return(NULL)
        }
    )

    if (is.null(resp)) {
        failures <- c(failures, acc)
        next
    }

    status <- httr2::resp_status(resp = resp)

    if (status != 200) {
        cat("HTTP", status, "\n")
        failures <- c(failures, acc)
        next
    }

    # Extract FASTA text
    fasta_text <- httr2::resp_body_string(resp = resp)

    # Parse: split into header and sequence
    lines <- strsplit(x = fasta_text, split = "\n", fixed = TRUE)[[1]]
    lines <- lines[nchar(lines) > 0]

    if (length(lines) < 2 || !startsWith(x = lines[1], prefix = ">")) {
        cat("INVALID FASTA\n")
        failures <- c(failures, acc)
        next
    }

    original_header <- lines[1]
    sequence <- paste(lines[-1], collapse = "")

    # Store with metadata
    sequences[[i]] <- list(
        organism = org,
        gene = gene,
        accession = acc,
        header = paste0(">", org, "_", gene),
        original_header = original_header,
        sequence = sequence
    )

    cat("OK (", nchar(sequence), " aa)\n", sep = "")

    # Rate limit delay with jitter
    delay <- REQUEST_DELAY_SEC + stats::runif(n = 1, min = -REQUEST_DELAY_JITTER, max = REQUEST_DELAY_JITTER)
    Sys.sleep(time = delay)
}

# Remove NULLs from failed requests
sequences <- Filter(f = Negate(is.null), x = sequences)

cat("\n=== FETCH SUMMARY ===\n")
cat("Successful:", length(sequences), "\n")
cat("Failed:", length(failures), "\n")
if (length(failures) > 0) {
    cat("Failed accessions:", paste(failures, collapse = ", "), "\n")
}

# ==============================================================================
# OUTPUT - WRITE PER-GENE FASTA FILES
# ==============================================================================

if (!dir.exists(OUTPUT_DIR)) {
    dir.create(path = OUTPUT_DIR, recursive = TRUE)
}

genes <- GENE_NAMES

cat("\n=== WRITING FASTA FILES ===\n")

for (g in genes) {
    gene_seqs <- Filter(f = function(x) x$gene == g, x = sequences)

    if (length(gene_seqs) == 0) {
        cat(g, ": no sequences, skipping\n")
        next
    }

    output_file <- file.path(OUTPUT_DIR, paste0(INPUT_PREFIX, g, ".fasta"))

    fasta_content <- character(length(gene_seqs) * 2)
    for (j in seq_along(gene_seqs)) {
        fasta_content[(j - 1) * 2 + 1] <- gene_seqs[[j]]$header
        fasta_content[(j - 1) * 2 + 2] <- gene_seqs[[j]]$sequence
    }

    writeLines(text = fasta_content, con = output_file)
    cat(g, ": wrote", length(gene_seqs), "sequences to", output_file, "\n")
}

cat("\nDone.\n")

# ==============================================================================
# PHASE 2: ALIGNMENT AND CONSERVATION
# ==============================================================================
# Adapted from stable_msa_conservation.R
# Aligns fetched sequences and calculates conservation relative to S. cerevisiae


# --- Phase 2 file paths (uses shared config from top) ---
if (!dir.exists(ALIGNMENT_OUTPUT_DIR)) {
    dir.create(path = ALIGNMENT_OUTPUT_DIR, recursive = TRUE)
}

align_input_files <- file.path(FASTA_OUTPUT_DIR, paste0(INPUT_PREFIX, GENE_NAMES, ".fasta"))
align_output_aligned <- file.path(ALIGNMENT_OUTPUT_DIR, paste0(INPUT_PREFIX, GENE_NAMES, "_aligned.fasta"))
align_output_conservation <- file.path(ALIGNMENT_OUTPUT_DIR, paste0(INPUT_PREFIX, GENE_NAMES, "_conservation.tsv"))

# Check inputs exist
align_missing <- align_input_files[!file.exists(align_input_files)]
if (length(align_missing) > 0) {
    message("Missing input files:\n", paste(align_missing, collapse = "\n"))
    stop("Fetch phase did not produce expected FASTA files.")
}

# Pre-allocate summary table
align_summary_df <- data.frame(
    gene = GENE_NAMES,
    n_sequences = integer(length(GENE_NAMES)),
    alignment_length = integer(length(GENE_NAMES)),
    scer_length = integer(length(GENE_NAMES)),
    mean_conservation = numeric(length(GENE_NAMES)),
    stringsAsFactors = FALSE
)

cat("\n=== PHASE 2: ALIGNMENT & CONSERVATION ===\n")
cat("Genes:", length(GENE_NAMES), "\n")
cat("Input:", FASTA_OUTPUT_DIR, "\n")
cat("Output:", ALIGNMENT_OUTPUT_DIR, "\n\n")

for (i in seq_along(GENE_NAMES)) {

    gene <- GENE_NAMES[i]
    cat("[", i, "/", length(GENE_NAMES), "] ", gene, "\n", sep = "")

    # --- Load sequences ---
    seqs <- Biostrings::readAAStringSet(filepath = align_input_files[i])
    n_seqs <- length(seqs)
    cat("  Loaded:", n_seqs, "sequences\n")

    # --- Align ---
    cat("  Aligning...")
    aligned <- DECIPHER::AlignSeqs(myXStringSet = seqs, verbose = FALSE)
    cat(" done\n")

    Biostrings::writeXStringSet(x = aligned, filepath = align_output_aligned[i])


    aln_widths <- unique(Biostrings::width(aligned))
    if (length(aln_widths) != 1) {
        stop("Alignment for ", gene, " has inconsistent sequence widths: ",
             paste(aln_widths, collapse = ", "))
    }
    aln_length <- aln_widths
    cat("  Alignment:", n_seqs, "x", aln_length, "\n")

    # --- Find reference sequence ---
    seq_names <- names(aligned)
    ref_idx <- grep(pattern = REFERENCE_PATTERN, x = seq_names)

    if (length(ref_idx) == 0) {
        cat("  WARNING: Reference not found, skipping conservation\n\n")
        align_summary_df$n_sequences[i] <- n_seqs
        align_summary_df$alignment_length[i] <- aln_length
        align_summary_df$scer_length[i] <- NA
        align_summary_df$mean_conservation[i] <- NA
        next
    }

    if (length(ref_idx) > 1) {
        cat("  WARNING: Multiple references, using first\n")
        ref_idx <- ref_idx[1]
    }

    # --- Conservation calculation ---
    aln_matrix <- as.matrix(aligned)
    ref_seq <- aln_matrix[ref_idx, ]

    ungapped_cols <- which(ref_seq != "-")
    scer_length <- length(ungapped_cols)

    conservation_df <- data.frame(
        position = seq_len(scer_length),
        residue = ref_seq[ungapped_cols],
        conservation_pct = numeric(scer_length),
        stringsAsFactors = FALSE
    )

    other_rows <- seq_len(n_seqs)[-ref_idx]
    other_matrix <- aln_matrix[other_rows, , drop = FALSE]

    for (j in seq_along(ungapped_cols)) {
        col_idx <- ungapped_cols[j]
        ref_residue <- ref_seq[col_idx]
        col_values <- other_matrix[, col_idx]

        non_gap <- col_values != "-"
        n_non_gap <- sum(non_gap)

        if (n_non_gap > 0) {
            n_identical <- sum(col_values[non_gap] == ref_residue)
            conservation_df$conservation_pct[j] <- (n_identical / n_non_gap) * 100
        } else {
            conservation_df$conservation_pct[j] <- NA
        }
    }

    utils::write.table(
        x = conservation_df,
        file = align_output_conservation[i],
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    mean_cons <- mean(conservation_df$conservation_pct, na.rm = TRUE)
    cat("  Conservation: mean =", round(mean_cons, 1), "% across", scer_length, "positions\n\n")

    align_summary_df$n_sequences[i] <- n_seqs
    align_summary_df$alignment_length[i] <- aln_length
    align_summary_df$scer_length[i] <- scer_length
    align_summary_df$mean_conservation[i] <- mean_cons
}

# --- Phase 2 Summary ---
align_summary_path <- file.path(ALIGNMENT_OUTPUT_DIR, paste0(INPUT_PREFIX, "summary.tsv"))
utils::write.table(
    x = align_summary_df,
    file = align_summary_path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)

cat("=== PHASE 2 COMPLETE ===\n")
cat("Summary saved to:", align_summary_path, "\n\n")
print(align_summary_df)

# ==============================================================================
# SESSION INFO
# ==============================================================================

session_info_path <- file.path(ALIGNMENT_OUTPUT_DIR, paste0(INPUT_PREFIX, "fetch_align_session_info.txt"))
session_info <- utils::capture.output(utils::sessionInfo())
writeLines(text = session_info, con = session_info_path)
cat("Session info:", session_info_path, "\n")
