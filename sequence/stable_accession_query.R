#!/usr/bin/env Rscript
# stable_accession_query.R
# Build and execute queries to fetch protein sequences from UniProt/NCBI

# ==============================================================================
# CONFIGURATION
# ==============================================================================

ROOT_DIRECTORY <- system("git rev-parse --show-toplevel", intern = TRUE)
INPUT_PATH <- file.path(ROOT_DIRECTORY, "sequence")

UNIPROT_BASE_URL <- "https://rest.uniprot.org/uniprotkb"
NCBI_EFETCH_URL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

INPUT_TSV <- file.path(INPUT_PATH, "stable_protein_accessions.tsv")
OUTPUT_DIR <- "~/data/protein_files/"
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

genes <- unique(accessions$gene)

cat("\n=== WRITING FASTA FILES ===\n")

for (g in genes) {
    gene_seqs <- Filter(f = function(x) x$gene == g, x = sequences)

    if (length(gene_seqs) == 0) {
        cat(g, ": no sequences, skipping\n")
        next
    }

    output_file <- file.path(OUTPUT_DIR, paste0("stable_", g, ".fasta"))

    fasta_content <- character(length(gene_seqs) * 2)
    for (j in seq_along(gene_seqs)) {
        fasta_content[(j - 1) * 2 + 1] <- gene_seqs[[j]]$header
        fasta_content[(j - 1) * 2 + 2] <- gene_seqs[[j]]$sequence
    }

    writeLines(text = fasta_content, con = output_file)
    cat(g, ": wrote", length(gene_seqs), "sequences to", output_file, "\n")
}

cat("\nDone.\n")
