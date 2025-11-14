# Date: 2025-11-14
# Notes:
#   Had to adjust plan to analyze genome with gff file.
#   Need to download those files from sgd.
#   Files were downloaded manually from SGD.
#   For now, use this script to analyze fasta files from SGD
#   After downloading to Downloads folder, run:
#     $for file in *.fsa;do mv $file "$HOME/data/fasta_files/" ; done
#   Switch mv to echo to perform dry-run.
# Purpose:
#   Find nggs in score in fasta files.
# Dependencies:
#   Assumes fasta files +/- 1000 flanking bp from sgd for gene of interest.
#   See libraries loaded for additional dependencies.
#   R 4.2
# ============================================================================
# CRISPR gRNA DESIGN - SGD FASTA FILES WITH CDS OVERLAP ANALYSIS
# Simplified pipeline for pre-downloaded gene sequences (ñ1000bp from SGD)
# ============================================================================

# ============================================================================
# REQUIRED LIBRARIES
# ============================================================================
library(Biostrings)
library(GenomicRanges)
library(IRanges)

# ============================================================================
# CONFIGURATION: Paths and parameters
# ============================================================================

# Input directory containing .fsa files from SGD
HOME_DIRECTORY <- Sys.getenv("HOME")
FASTA_FILE_DIRECTORY_NAME <- "data/fasta_files"
INPUT_DIR_PATH_chr <- file.path(HOME_DIRECTORY, FASTA_FILE_DIRECTORY_NAME)

# Output file

OUTPUT_FILENAME_chr <- "grna_results.tsv"
OUTPUT_FILE_PATH_chr <- file.path(
  HOME_DIRECTORY, FASTA_FILE_DIRECTORY_NAME,
  OUTPUT_FILENAME_chr
)

# CDS position constants (SGD W303 format)
CDS_START_EXPECTED_pos_int <- 1001  # Start codon always at position 1001
CDS_END_OFFSET_bp_int <- 1000       # Stop codon always 1001bp from end

# Search region extension around CDS
UPSTREAM_EXTENSION_bp_int <- 100
DOWNSTREAM_EXTENSION_bp_int <- 100

# Guide RNA parameters
GUIDE_LENGTH_bp_int <- 20
PAM_SEQUENCE_chr <- "NGG"

# ============================================================================
# CONFIGURATION: Genetic code table
# ============================================================================

GENETIC_CODE_nls <- list(
  TTT = "F", TTC = "F", TTA = "L", TTG = "L",
  TCT = "S", TCC = "S", TCA = "S", TCG = "S",
  TAT = "Y", TAC = "Y", TAA = "*", TAG = "*",
  TGT = "C", TGC = "C", TGA = "*", TGG = "W",
  CTT = "L", CTC = "L", CTA = "L", CTG = "L",
  CCT = "P", CCC = "P", CCA = "P", CCG = "P",
  CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
  CGT = "R", CGC = "R", CGA = "R", CGG = "R",
  ATT = "I", ATC = "I", ATA = "I", ATG = "M",
  ACT = "T", ACC = "T", ACA = "T", ACG = "T",
  AAT = "N", AAC = "N", AAA = "K", AAG = "K",
  AGT = "S", AGC = "S", AGA = "R", AGG = "R",
  GTT = "V", GTC = "V", GTA = "V", GTG = "V",
  GCT = "A", GCC = "A", GCA = "A", GCG = "A",
  GAT = "D", GAC = "D", GAA = "E", GAG = "E",
  GGT = "G", GGC = "G", GGA = "G", GGG = "G"
)

# ============================================================================
# STEP 1: List all .fsa files in input directory
# ============================================================================

cat("Scanning directory for .fsa files:", INPUT_DIR_PATH_chr, "\n")

fasta_files_chr <- list.files(path = INPUT_DIR_PATH_chr,
                               pattern = "\\.fsa$",
                               full.names = TRUE)

cat("Found", length(fasta_files_chr), "FASTA files\n")
cat("Files:", paste(basename(fasta_files_chr), collapse = ", "), "\n\n")

# ============================================================================
# STEP 2: Test reading first file and verify CDS boundaries
# ============================================================================

cat("Testing first file:", basename(fasta_files_chr[1]), "\n")

# Read sequence
test_seq_dss <- Biostrings::readDNAStringSet(filepath = fasta_files_chr[1],
                                              format = "fasta")

# Extract gene name from header (first word)
header_chr <- names(test_seq_dss)[1]
gene_name_chr <- strsplit(x = header_chr, split = " ")[[1]][1]

cat("Gene name:", gene_name_chr, "\n")

# Get sequence
seq_dna <- test_seq_dss[[1]]
seq_length_bp_int <- length(seq_dna)

cat("Sequence length:", seq_length_bp_int, "bp\n")

# ============================================================================
# STEP 3: Verify start codon at position 1001
# ============================================================================

start_codon_chr <- as.character(Biostrings::subseq(x = seq_dna,
                                                     start = CDS_START_EXPECTED_pos_int,
                                                     end = CDS_START_EXPECTED_pos_int + 2))

cat("Start codon at position", CDS_START_EXPECTED_pos_int, ":", start_codon_chr, "\n")

if (start_codon_chr != "ATG") {
  stop("ERROR: Start codon is not ATG at position ", CDS_START_EXPECTED_pos_int,
       ". Found: ", start_codon_chr)
}

cat(" Start codon verified: ATG\n")

# ============================================================================
# STEP 4: Verify stop codon at expected position
# ============================================================================

# Calculate stop codon position
cds_end_pos_int <- seq_length_bp_int - CDS_END_OFFSET_bp_int
stop_codon_start_pos_int <- cds_end_pos_int - 2

stop_codon_chr <- as.character(Biostrings::subseq(x = seq_dna,
                                                    start = stop_codon_start_pos_int,
                                                    end = cds_end_pos_int))

cat("Stop codon at position", stop_codon_start_pos_int, "-", cds_end_pos_int, ":", 
    stop_codon_chr, "\n")

valid_stop_codons_chr <- c("TAA", "TAG", "TGA")

if (!stop_codon_chr %in% valid_stop_codons_chr) {
  stop("ERROR: Stop codon is not TAA/TAG/TGA at expected position. Found: ", 
       stop_codon_chr)
}

cat(" Stop codon verified:", stop_codon_chr, "\n")

# ============================================================================
# STEP 5: Calculate and verify CDS length
# ============================================================================

cds_length_bp_int <- cds_end_pos_int - CDS_START_EXPECTED_pos_int + 1

cat("CDS length:", cds_length_bp_int, "bp (", cds_length_bp_int / 3, "codons)\n")

# Verify CDS length is multiple of 3
if (cds_length_bp_int %% 3 != 0) {
  stop("ERROR: CDS length is not a multiple of 3: ", cds_length_bp_int)
}

cat(" CDS length is valid (multiple of 3)\n\n")

cat("=== VERIFICATION COMPLETE ===\n")
cat("Gene:", gene_name_chr, "\n")
cat("CDS: position", CDS_START_EXPECTED_pos_int, "to", cds_end_pos_int, 
    "(", cds_length_bp_int, "bp)\n")
