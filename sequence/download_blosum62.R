#!/usr/bin/env Rscript
# download_blosum62.R
# Download and cache BLOSUM62 matrix from NCBI with provenance metadata
#
# Usage: Rscript download_blosum62.R [--force]
#   --force: Re-download even if cached file exists

# ==============================================================================
# CONFIGURATION
# ==============================================================================

BLOSUM62_URL <- "https://ftp.ncbi.nih.gov/blast/matrices/BLOSUM62"
OUTPUT_DIR <- "~/data/protein_files"
BLOSUM62_FILE <- file.path(OUTPUT_DIR, "BLOSUM62.txt")
BLOSUM62_RDS <- file.path(OUTPUT_DIR, "BLOSUM62.rds")
METADATA_FILE <- file.path(OUTPUT_DIR, "BLOSUM62_metadata.txt")

# Expected properties for validation
EXPECTED_NROW <- 24
EXPECTED_NCOL <- 24
EXPECTED_AA <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","B","Z","X","*")

# ==============================================================================
# PARSE ARGUMENTS
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
force_download <- "--force" %in% args

# ==============================================================================
# SETUP
# ==============================================================================

OUTPUT_DIR <- path.expand(OUTPUT_DIR)
BLOSUM62_FILE <- path.expand(BLOSUM62_FILE)
BLOSUM62_RDS <- path.expand(BLOSUM62_RDS)
METADATA_FILE <- path.expand(METADATA_FILE)

if (!dir.exists(OUTPUT_DIR)) {
    dir.create(path = OUTPUT_DIR, recursive = TRUE)
}

cat("=== BLOSUM62 DOWNLOAD ===\n")
cat("Source:", BLOSUM62_URL, "\n")
cat("Output:", BLOSUM62_FILE, "\n\n")

# ==============================================================================
# CHECK CACHE
# ==============================================================================

if (file.exists(BLOSUM62_RDS) && file.exists(METADATA_FILE) && !force_download) {
    cat("Cached files found. Use --force to re-download.\n")
    cat("Loading cached matrix to verify...\n")
    
    BLOSUM62 <- readRDS(file = BLOSUM62_RDS)
    cat("  Dimensions:", nrow(BLOSUM62), "x", ncol(BLOSUM62), "\n")
    cat("  Row names:", paste(head(rownames(BLOSUM62), 5), collapse = ", "), ", ...\n")
    cat("\nCache valid. Exiting.\n")
    quit(save = "no", status = 0)
}

# ==============================================================================
# DOWNLOAD
# ==============================================================================

cat("Downloading BLOSUM62 matrix...\n")
download_time <- Sys.time()

# Download to temporary file first
temp_file <- tempfile(fileext = ".txt")
download_result <- utils::download.file(
    url = BLOSUM62_URL,
    destfile = temp_file,
    method = "auto",
    quiet = FALSE
)

if (download_result != 0) {
    stop("Download failed with code: ", download_result)
}

cat("Download complete.\n\n")

# ==============================================================================
# COMPUTE MD5 OF RAW FILE
# ==============================================================================

raw_md5 <- tools::md5sum(temp_file)
names(raw_md5) <- NULL
cat("Raw file MD5:", raw_md5, "\n\n")

# ==============================================================================
# PARSE MATRIX
# ==============================================================================

cat("Parsing matrix...\n")

# Read all lines
raw_lines <- readLines(con = temp_file)

# Filter out comment lines (start with #) and empty lines
data_lines <- raw_lines[!grepl(pattern = "^#", x = raw_lines) & nchar(trimws(raw_lines)) > 0]

# First data line contains column headers
header_line <- data_lines[1]
col_names <- strsplit(x = trimws(header_line), split = "\\s+")[[1]]

# Remaining lines contain row data
data_lines <- data_lines[-1]

# Parse each row
n_rows <- length(data_lines)
matrix_data <- matrix(data = NA_integer_, nrow = n_rows, ncol = length(col_names))
row_names <- character(n_rows)

for (i in seq_along(data_lines)) {
    fields <- strsplit(x = trimws(data_lines[i]), split = "\\s+")[[1]]
    row_names[i] <- fields[1]
    matrix_data[i, ] <- as.integer(fields[-1])
}

# Assemble matrix with dimnames
BLOSUM62 <- matrix_data
rownames(BLOSUM62) <- row_names
colnames(BLOSUM62) <- col_names

cat("Parsed matrix:", nrow(BLOSUM62), "x", ncol(BLOSUM62), "\n\n")

# ==============================================================================
# VALIDATION
# ==============================================================================

cat("Validating matrix...\n")

# Check dimensions
stopifnot(
    "Unexpected row count" = nrow(BLOSUM62) == EXPECTED_NROW,
    "Unexpected column count" = ncol(BLOSUM62) == EXPECTED_NCOL
)
cat("  Dimensions: OK (", EXPECTED_NROW, "x", EXPECTED_NCOL, ")\n", sep = "")

# Check row/column names match expected amino acids
stopifnot(
    "Row names mismatch" = all(rownames(BLOSUM62) == EXPECTED_AA),
    "Column names mismatch" = all(colnames(BLOSUM62) == EXPECTED_AA)
)
cat("  Amino acid labels: OK\n")

# Check symmetry (BLOSUM62 should be symmetric)
stopifnot(
    "Matrix not symmetric" = all(BLOSUM62 == t(BLOSUM62))
)
cat("  Symmetry: OK\n")

# Check diagonal values are positive (self-similarity)
diag_values <- diag(BLOSUM62)
stopifnot(
    "Diagonal contains non-positive values for standard AAs" = all(diag_values[1:20] > 0)
)
cat("  Diagonal (standard AAs): OK (all positive)\n")

# Check known values (spot check)
stopifnot(
    "A-A should be 4" = BLOSUM62["A", "A"] == 4,
    "W-W should be 11" = BLOSUM62["W", "W"] == 11,
    "W-C should be -2" = BLOSUM62["W", "C"] == -2
)
cat("  Spot check values: OK\n")

cat("\nAll validations passed.\n\n")

# ==============================================================================
# SAVE OUTPUTS
# ==============================================================================

# Copy raw file to permanent location
file.copy(from = temp_file, to = BLOSUM62_FILE, overwrite = TRUE)
cat("Saved raw text:", BLOSUM62_FILE, "\n")

# Save as RDS for fast loading
saveRDS(object = BLOSUM62, file = BLOSUM62_RDS)
cat("Saved RDS:", BLOSUM62_RDS, "\n")

# ==============================================================================
# WRITE METADATA
# ==============================================================================

metadata_lines <- c(
    "# BLOSUM62 Matrix Provenance",
    paste0("source_url: ", BLOSUM62_URL),
    paste0("download_date: ", format(download_time, "%Y-%m-%d %H:%M:%S %Z")),
    paste0("raw_file_md5: ", raw_md5),
    paste0("dimensions: ", nrow(BLOSUM62), "x", ncol(BLOSUM62)),
    paste0("row_names: ", paste(rownames(BLOSUM62), collapse = ",")),
    paste0("col_names: ", paste(colnames(BLOSUM62), collapse = ",")),
    paste0("validation_passed: TRUE"),
    paste0("r_version: ", R.version.string)
)

writeLines(text = metadata_lines, con = METADATA_FILE)
cat("Saved metadata:", METADATA_FILE, "\n")

# ==============================================================================
# CLEANUP
# ==============================================================================

unlink(temp_file)

cat("\n=== COMPLETE ===\n")
cat("Matrix ready for use. Load with:\n")
cat('  BLOSUM62 <- readRDS("', BLOSUM62_RDS, '")\n', sep = "")
