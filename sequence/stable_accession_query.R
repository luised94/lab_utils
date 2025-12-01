#!/usr/bin/env Rscript
# stable_accession_query.R
# Build and optionally execute queries to fetch protein sequences from UniProt/NCBI

# ==============================================================================
# CONFIGURATION
# ==============================================================================
ROOT_DIRECTORY <- system("git rev-parse --show-toplevel", intern = TRUE)
INPUT_PATH <- file.path(ROOT_DIRECTORY, "sequence")

UNIPROT_BASE_URL <- "https://rest.uniprot.org/uniprotkb"
NCBI_EFETCH_URL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

INPUT_TSV <- file.path(INPUT_PATH, "stable_protein_accessions.tsv")
OUTPUT_DIR <- "protein_files"
QUERIES_OUTPUT <- "stable_fetch_queries.tsv"

# Request settings
REQUEST_DELAY_SEC <- 0.35
REQUEST_TIMEOUT_SEC <- 30

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

# Print all URLs grouped by gene (for MSA preparation)
cat("=== ALL URLS BY GENE ===\n\n")

genes <- unique(accessions$gene)

for (g in genes) {
    gene_rows <- accessions[accessions$gene == g, ]
    cat("---", g, "(n =", nrow(gene_rows), ") ---\n")
    for (j in seq_len(nrow(gene_rows))) {
        cat(gene_rows$organism[j], "\t", gene_rows$url[j], "\n")
    }
    cat("\n")
}
