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

# ============================================================================
# STEP 6: Calculate search region (CDS ñ 100bp)
# ============================================================================

search_start_pos_int <- max(1, CDS_START_EXPECTED_pos_int - UPSTREAM_EXTENSION_bp_int)
search_end_pos_int <- min(seq_length_bp_int, 
                           cds_end_pos_int + DOWNSTREAM_EXTENSION_bp_int)

cat("Search region: position", search_start_pos_int, "to", search_end_pos_int,
    "(", search_end_pos_int - search_start_pos_int + 1, "bp)\n\n")

# Extract search sequence
search_seq_dna <- Biostrings::subseq(x = seq_dna,
                                      start = search_start_pos_int,
                                      end = search_end_pos_int)

cat("Search sequence length:", length(search_seq_dna), "bp\n")

# ============================================================================
# STEP 7: Find PAM sites on forward strand
# ============================================================================

cat("\nSearching for PAM sites (", PAM_SEQUENCE_chr, ") on forward strand...\n", sep = "")

pam_matches_fwd_vws <- Biostrings::matchPattern(pattern = PAM_SEQUENCE_chr,
                                                  subject = search_seq_dna,
                                                  fixed = FALSE)

cat("Found", length(pam_matches_fwd_vws), "PAM sites on forward strand\n")

# ============================================================================
# STEP 8: Find PAM sites on reverse strand
# ============================================================================

cat("Searching for PAM sites on reverse strand...\n")

# Get reverse complement of search sequence
search_seq_rc_dna <- Biostrings::reverseComplement(x = search_seq_dna)

pam_matches_rev_vws <- Biostrings::matchPattern(pattern = PAM_SEQUENCE_chr,
                                                  subject = search_seq_rc_dna,
                                                  fixed = FALSE)

cat("Found", length(pam_matches_rev_vws), "PAM sites on reverse strand\n")

# ============================================================================
# STEP 9: Extract guide sequences for forward strand PAMs
# ============================================================================

cat("\nExtracting guide sequences for forward strand...\n")

guides_fwd_lst <- list()
search_seq_length_bp_int <- length(search_seq_dna)

for (i in seq_along(pam_matches_fwd_vws)) {
  # Get PAM position in search sequence
  pam_start_search_int <- BiocGenerics::start(pam_matches_fwd_vws)[i]
  pam_end_search_int <- BiocGenerics::end(pam_matches_fwd_vws)[i]

  # Calculate guide position (20bp upstream of PAM)
  guide_start_search_int <- pam_start_search_int - GUIDE_LENGTH_bp_int
  guide_end_search_int <- pam_start_search_int - 1

  # Check if we have enough sequence for full-length guide
  if (guide_start_search_int >= 1) {
    # Extract sequences
    guide_seq_chr <- as.character(Biostrings::subseq(x = search_seq_dna,
                                                       start = guide_start_search_int,
                                                       end = guide_end_search_int))

    pam_seq_chr <- as.character(Biostrings::subseq(x = search_seq_dna,
                                                     start = pam_start_search_int,
                                                     end = pam_end_search_int))

    full_seq_chr <- as.character(Biostrings::subseq(x = search_seq_dna,
                                                      start = guide_start_search_int,
                                                      end = pam_end_search_int))

    # Convert to original sequence coordinates
    guide_start_orig_int <- search_start_pos_int + guide_start_search_int - 1
    guide_end_orig_int <- search_start_pos_int + guide_end_search_int - 1
    pam_start_orig_int <- search_start_pos_int + pam_start_search_int - 1
    pam_end_orig_int <- search_start_pos_int + pam_end_search_int - 1

    # Store guide information
    guides_fwd_lst[[length(guides_fwd_lst) + 1]] <- list(
      gene_name = gene_name_chr,
      guide_strand = "+",
      guide_sequence = guide_seq_chr,
      pam_sequence = pam_seq_chr,
      full_sequence = full_seq_chr,
      guide_start_pos = guide_start_orig_int,
      guide_end_pos = guide_end_orig_int,
      pam_start_pos = pam_start_orig_int,
      pam_end_pos = pam_end_orig_int
    )
  }
}

cat("Extracted", length(guides_fwd_lst), "full-length guides from forward strand\n")

# ============================================================================
# STEP 10: Extract guide sequences for reverse strand PAMs
# ============================================================================

cat("Extracting guide sequences for reverse strand...\n")

guides_rev_lst <- list()

for (i in seq_along(pam_matches_rev_vws)) {
  # Get PAM position in reverse complement
  pam_start_rc_int <- BiocGenerics::start(pam_matches_rev_vws)[i]
  pam_end_rc_int <- BiocGenerics::end(pam_matches_rev_vws)[i]

  # Calculate guide position in reverse complement
  guide_start_rc_int <- pam_start_rc_int - GUIDE_LENGTH_bp_int
  guide_end_rc_int <- pam_start_rc_int - 1

  # Check if we have enough sequence for full-length guide
  if (guide_start_rc_int >= 1) {
    # Extract from reverse complement
    guide_seq_rc_chr <- as.character(Biostrings::subseq(x = search_seq_rc_dna,
                                                          start = guide_start_rc_int,
                                                          end = guide_end_rc_int))

    pam_seq_rc_chr <- as.character(Biostrings::subseq(x = search_seq_rc_dna,
                                                        start = pam_start_rc_int,
                                                        end = pam_end_rc_int))

    full_seq_rc_chr <- as.character(Biostrings::subseq(x = search_seq_rc_dna,
                                                         start = guide_start_rc_int,
                                                         end = pam_end_rc_int))


    # Convert RC coordinates to original sequence coordinates
    pam_end_orig_int <- search_start_pos_int + search_seq_length_bp_int - pam_start_rc_int
    pam_start_orig_int <- search_start_pos_int + search_seq_length_bp_int - pam_end_rc_int
    guide_end_orig_int <- search_start_pos_int + search_seq_length_bp_int - guide_start_rc_int
    guide_start_orig_int <- search_start_pos_int + search_seq_length_bp_int - guide_end_rc_int

    # Store guide information
    guides_rev_lst[[length(guides_rev_lst) + 1]] <- list(
      gene_name = gene_name_chr,
      guide_strand = "-",
      guide_sequence = guide_seq_rc_chr,
      pam_sequence = pam_seq_rc_chr,
      full_sequence = full_seq_rc_chr,
      guide_start_pos = guide_start_orig_int,
      guide_end_pos = guide_end_orig_int,
      pam_start_pos = pam_start_orig_int,
      pam_end_pos = pam_end_orig_int
    )
  }
}

cat("Extracted", length(guides_rev_lst), "full-length guides from reverse strand\n")

# ============================================================================
# STEP 11: Combine guides from both strands
# ============================================================================

all_guides_lst <- c(guides_fwd_lst, guides_rev_lst)

cat("\n=== GUIDE EXTRACTION COMPLETE ===\n")
cat("Total guides found:", length(all_guides_lst), "\n")
cat("Forward strand:", length(guides_fwd_lst), "\n")
cat("Reverse strand:", length(guides_rev_lst), "\n")
