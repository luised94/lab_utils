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

if (length(fasta_files_chr) == 0){
  stop("No fasta files found in INPUT_DIR_PATH_chr: ", INPUT_DIR_PATH_chr)
}

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

    # === VALIDATION ASSERTIONS ===
    if (nchar(guide_seq_chr) != 20) {
      stop("ERROR: Guide sequence length is ", nchar(guide_seq_chr), 
           " (expected 20) at position ", guide_start_orig_int)
    }

    if (nchar(pam_seq_chr) != 3) {
      stop("ERROR: PAM sequence length is ", nchar(pam_seq_chr), 
           " (expected 3) at position ", pam_start_orig_int)
    }

    if (substr(pam_seq_chr, 2, 3) != "GG") {
      stop("ERROR: PAM sequence does not end with GG. Found: ", pam_seq_chr,
           " at position ", pam_start_orig_int)
    }

    if (nchar(full_seq_chr) != 23) {
      stop("ERROR: Full sequence length is ", nchar(full_seq_chr),
           " (expected 23) at position ", guide_start_orig_int)
    }

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


    # === VALIDATION ASSERTIONS ===
    if (nchar(guide_seq_rc_chr) != 20) {
      stop("ERROR: Guide sequence length is ", nchar(guide_seq_rc_chr),
           " (expected 20) for reverse strand at RC position ", guide_start_rc_int)
    }

    if (nchar(pam_seq_rc_chr) != 3) {
      stop("ERROR: PAM sequence length is ", nchar(pam_seq_rc_chr),
           " (expected 3) for reverse strand at RC position ", pam_start_rc_int)
    }

    if (substr(pam_seq_rc_chr, 2, 3) != "GG") {
      stop("ERROR: PAM sequence does not end with GG. Found: ", pam_seq_rc_chr,
           " for reverse strand at RC position ", pam_start_rc_int)
    }

    if (nchar(full_seq_rc_chr) != 23) {
      stop("ERROR: Full sequence length is ", nchar(full_seq_rc_chr),
           " (expected 23) for reverse strand at RC position ", guide_start_rc_int)
    }

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

# ============================================================================
# STEP 12: Analyze CDS overlap and synonymous mutations for each guide
# ============================================================================

cat("\nAnalyzing CDS overlap and synonymous mutations...\n")

# CDS boundaries (already calculated)
cds_start_pos_int <- CDS_START_EXPECTED_pos_int  # 1001
cds_end_pos_int <- cds_end_pos_int  # Already calculated
cds_phase_int <- 0  # Phase is always 0 for these SGD files (starts at ATG)
# ============================================================================
# CDS overlap and synonymous mutation analysis
# ============================================================================

# Extract CDS sequence once
cds_seq_dna <- Biostrings::subseq(x = seq_dna,
                                   start = cds_start_pos_int,
                                   end = cds_end_pos_int)
cds_seq_chr <- as.character(cds_seq_dna)

cat("\nAnalyzing CDS overlap and synonymous mutations...\n")

guides_analyzed_int <- 0

for (guide_idx in seq_along(all_guides_lst)) {
  guide_info_lst <- all_guides_lst[[guide_idx]]

  # Check overlap
  guide_region_start_int <- guide_info_lst$guide_start_pos
  guide_region_end_int <- guide_info_lst$pam_end_pos

  overlaps_cds_lgl <- (guide_region_start_int <= cds_end_pos_int && 
                       guide_region_end_int >= cds_start_pos_int)

  guide_info_lst$overlaps_cds <- overlaps_cds_lgl

  # If it overlaps cds, must determine if pam can be mutated easily.
  if (overlaps_cds_lgl) {
    guide_info_lst$cds_start <- cds_start_pos_int
    guide_info_lst$cds_end <- cds_end_pos_int

    # =======================================================================
    # Determine which PAM positions to check based on strand
    # =======================================================================

    if (guide_info_lst$guide_strand == "+") {
      # Forward strand: PAM is NGG, check positions 2-3 (the GG)
      pam_check_positions_int <- 2:3
      pam_for_analysis_chr <- guide_info_lst$pam_sequence

    } else {
      # Reverse strand: stored PAM is RC form
      # Need to analyze forward strand, check positions 1-2
      pam_check_positions_int <- 1:2

      # Reverse complement to get forward strand representation
      pam_for_analysis_chr <- as.character(
        Biostrings::reverseComplement(
          Biostrings::DNAString(guide_info_lst$pam_sequence)
        )
      )
    }

    # =======================================================================
    # Get PAM positions and check if can mutate synonymously
    # =======================================================================

    pam_positions_int <- guide_info_lst$pam_start_pos:guide_info_lst$pam_end_pos
    pam_can_mutate_lgl <- FALSE

    for (pam_idx in pam_check_positions_int) {
      pam_abs_pos_int <- pam_positions_int[pam_idx]

      # Position in CDS (1-based)
      pos_in_cds_int <- pam_abs_pos_int - cds_start_pos_int + 1

      # Which codon and position within codon
      codon_number_int <- ceiling(pos_in_cds_int / 3)
      position_in_codon_int <- ((pos_in_cds_int - 1) %% 3) + 1

      # Extract codon
      codon_start_int <- (codon_number_int - 1) * 3 + 1
      codon_end_int <- codon_start_int + 2

      if (codon_end_int <= nchar(cds_seq_chr)) {
        codon_chr <- substr(cds_seq_chr, codon_start_int, codon_end_int)
        current_aa_chr <- GENETIC_CODE_nls[[codon_chr]]

        if (!is.null(current_aa_chr) && !is.na(current_aa_chr)) {
          current_nt_chr <- substr(codon_chr, position_in_codon_int, position_in_codon_int)

          # Test all mutations
          for (new_nt_chr in c("A", "T", "G", "C")) {
            if (new_nt_chr != current_nt_chr) {
              mutated_codon_chr <- codon_chr
              substr(mutated_codon_chr, position_in_codon_int, position_in_codon_int) <- new_nt_chr

              mutated_aa_chr <- GENETIC_CODE_nls[[mutated_codon_chr]]

              if (!is.null(mutated_aa_chr) && !is.na(mutated_aa_chr) && 
                  mutated_aa_chr == current_aa_chr) {
                pam_can_mutate_lgl <- TRUE
                break
              }
            }
          }
        }
      }

      if (pam_can_mutate_lgl) break  # Found one, no need to check other position
    }

    guide_info_lst$pam_gg_can_mutate <- pam_can_mutate_lgl

    # =======================================================================
    # Count wobble positions in guide
    # =======================================================================

    guide_positions_int <- guide_info_lst$guide_start_pos:guide_info_lst$guide_end_pos
    guide_in_cds_int <- guide_positions_int[
      guide_positions_int >= cds_start_pos_int & 
      guide_positions_int <= cds_end_pos_int
    ]

    if (length(guide_in_cds_int) > 0) {
      pos_in_cds_int <- guide_in_cds_int - cds_start_pos_int + 1
      position_in_codon_int <- ((pos_in_cds_int - 1) %% 3) + 1
      wobble_count_int <- sum(position_in_codon_int == 3)
    } else {
      wobble_count_int <- 0
    }

    guide_info_lst$guide_wobble_sites <- wobble_count_int

    # =======================================================================
    # Silencing difficulty
    # =======================================================================

    if (guide_info_lst$pam_gg_can_mutate) {
      guide_info_lst$silencing_difficulty <- "easy"
      guide_info_lst$silencing_strategy <- "PAM_mutation"
    } else if (wobble_count_int >= 3) {
      guide_info_lst$silencing_difficulty <- "moderate"
      guide_info_lst$silencing_strategy <- "guide_mutations"
    } else if (wobble_count_int >= 1) {
      guide_info_lst$silencing_difficulty <- "difficult"
      guide_info_lst$silencing_strategy <- "limited_guide_mutations"
    } else {
      guide_info_lst$silencing_difficulty <- "very_difficult"
      guide_info_lst$silencing_strategy <- "non_synonymous_required"
    }

  } else {
    # Non-coding
    guide_info_lst$cds_start <- NA
    guide_info_lst$cds_end <- NA
    guide_info_lst$pam_gg_can_mutate <- TRUE
    guide_info_lst$guide_wobble_sites <- NA
    guide_info_lst$silencing_difficulty <- "easy"
    guide_info_lst$silencing_strategy <- "non_coding_region"
  }

  all_guides_lst[[guide_idx]] <- guide_info_lst
  guides_analyzed_int <- guides_analyzed_int + 1
}

cat("Analyzed", guides_analyzed_int, "guides\n")

# Summary
cat("\n=== ANALYSIS SUMMARY ===\n")
cat("Guides overlapping CDS:", sum(sapply(all_guides_lst, function(x) x$overlaps_cds)), "\n")
cat("Difficulty distribution:\n")
print(table(sapply(all_guides_lst, function(x) x$silencing_difficulty)))
# ============================================================================
# STEP 14: Process all FASTA files
# ============================================================================

cat("\n=== PROCESSING ALL FILES ===\n")

# Storage for all genes
all_genes_results_lst <- list()

for (file_idx in seq_along(fasta_files_chr)) {
  fasta_file_path_chr <- fasta_files_chr[file_idx]

  cat("\nProcessing file", file_idx, "of", length(fasta_files_chr), ":",
      basename(fasta_file_path_chr), "\n")

  # Read sequence
  seq_dss <- Biostrings::readDNAStringSet(filepath = fasta_file_path_chr,
                                          format = "fasta")

  # Extract gene name
  header_chr <- names(seq_dss)[1]
  gene_name_chr <- strsplit(x = header_chr, split = " ")[[1]][1]

  cat("  Gene:", gene_name_chr, "\n")

  # Get sequence
  seq_dna <- seq_dss[[1]]
  seq_length_bp_int <- length(seq_dna)

  # ==========================================================================
  # Verify CDS boundaries
  # ==========================================================================

  # Check start codon
  start_codon_chr <- as.character(Biostrings::subseq(x = seq_dna,
                                                       start = CDS_START_EXPECTED_pos_int,
                                                       end = CDS_START_EXPECTED_pos_int + 2))

  if (start_codon_chr != "ATG") {
    warning(paste("Gene", gene_name_chr, ": Start codon is not ATG at position",
                  CDS_START_EXPECTED_pos_int, ". Found:", start_codon_chr, "- SKIPPING"))
    next
  }

  # Calculate stop codon position
  cds_end_pos_int <- seq_length_bp_int - CDS_END_OFFSET_bp_int
  stop_codon_start_pos_int <- cds_end_pos_int - 2

  stop_codon_chr <- as.character(Biostrings::subseq(x = seq_dna,
                                                      start = stop_codon_start_pos_int,
                                                      end = cds_end_pos_int))

  valid_stop_codons_chr <- c("TAA", "TAG", "TGA")

  if (!stop_codon_chr %in% valid_stop_codons_chr) {
    warning(paste("Gene", gene_name_chr, ": Stop codon is not valid at expected position.",
                  "Found:", stop_codon_chr, "- SKIPPING"))
    next
  }

  # Calculate CDS length
  cds_length_bp_int <- cds_end_pos_int - CDS_START_EXPECTED_pos_int + 1

  if (cds_length_bp_int %% 3 != 0) {
    warning(paste("Gene", gene_name_chr, ": CDS length is not a multiple of 3:",
                  cds_length_bp_int, "- SKIPPING"))
    next
  }

  cat("  CDS:", CDS_START_EXPECTED_pos_int, "-", cds_end_pos_int,
      "(", cds_length_bp_int, "bp,", cds_length_bp_int / 3, "codons)\n")

  # ==========================================================================
  # Calculate search region
  # ==========================================================================

  search_start_pos_int <- max(1, CDS_START_EXPECTED_pos_int - UPSTREAM_EXTENSION_bp_int)
  search_end_pos_int <- min(seq_length_bp_int,
                            cds_end_pos_int + DOWNSTREAM_EXTENSION_bp_int)

  search_seq_dna <- Biostrings::subseq(x = seq_dna,
                                        start = search_start_pos_int,
                                        end = search_end_pos_int)

  search_seq_length_bp_int <- length(search_seq_dna)

  # ==========================================================================
  # Find PAMs on both strands
  # ==========================================================================

  # Forward strand
  pam_matches_fwd_vws <- Biostrings::matchPattern(pattern = PAM_SEQUENCE_chr,
                                                    subject = search_seq_dna,
                                                    fixed = FALSE)

  # Reverse strand
  search_seq_rc_dna <- Biostrings::reverseComplement(x = search_seq_dna)
  pam_matches_rev_vws <- Biostrings::matchPattern(pattern = PAM_SEQUENCE_chr,
                                                    subject = search_seq_rc_dna,
                                                    fixed = FALSE)

  cat("  PAMs found: forward =", length(pam_matches_fwd_vws),
      ", reverse =", length(pam_matches_rev_vws), "\n")

  # ==========================================================================
  # Extract guides from forward strand
  # ==========================================================================

  guides_fwd_lst <- list()

  for (i in seq_along(pam_matches_fwd_vws)) {
    pam_start_search_int <- BiocGenerics::start(pam_matches_fwd_vws)[i]
    pam_end_search_int <- BiocGenerics::end(pam_matches_fwd_vws)[i]

    guide_start_search_int <- pam_start_search_int - GUIDE_LENGTH_bp_int
    guide_end_search_int <- pam_start_search_int - 1

    if (guide_start_search_int >= 1) {
      guide_seq_chr <- as.character(Biostrings::subseq(x = search_seq_dna,
                                                         start = guide_start_search_int,
                                                         end = guide_end_search_int))

      pam_seq_chr <- as.character(Biostrings::subseq(x = search_seq_dna,
                                                       start = pam_start_search_int,
                                                       end = pam_end_search_int))

      full_seq_chr <- as.character(Biostrings::subseq(x = search_seq_dna,
                                                        start = guide_start_search_int,
                                                        end = pam_end_search_int))

      # Validate
      if (nchar(guide_seq_chr) != 20 || nchar(pam_seq_chr) != 3 || 
          substr(pam_seq_chr, 2, 3) != "GG" || nchar(full_seq_chr) != 23) {
        next
      }

      # Convert to original coordinates
      guide_start_orig_int <- search_start_pos_int + guide_start_search_int - 1
      guide_end_orig_int <- search_start_pos_int + guide_end_search_int - 1
      pam_start_orig_int <- search_start_pos_int + pam_start_search_int - 1
      pam_end_orig_int <- search_start_pos_int + pam_end_search_int - 1

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

  # ==========================================================================
  # Extract guides from reverse strand
  # ==========================================================================

  guides_rev_lst <- list()

  for (i in seq_along(pam_matches_rev_vws)) {
    pam_start_rc_int <- BiocGenerics::start(pam_matches_rev_vws)[i]
    pam_end_rc_int <- BiocGenerics::end(pam_matches_rev_vws)[i]

    guide_start_rc_int <- pam_start_rc_int - GUIDE_LENGTH_bp_int
    guide_end_rc_int <- pam_start_rc_int - 1

    if (guide_start_rc_int >= 1) {
      guide_seq_rc_chr <- as.character(Biostrings::subseq(x = search_seq_rc_dna,
                                                            start = guide_start_rc_int,
                                                            end = guide_end_rc_int))

      pam_seq_rc_chr <- as.character(Biostrings::subseq(x = search_seq_rc_dna,
                                                          start = pam_start_rc_int,
                                                          end = pam_end_rc_int))

      full_seq_rc_chr <- as.character(Biostrings::subseq(x = search_seq_rc_dna,
                                                           start = guide_start_rc_int,
                                                           end = pam_end_rc_int))

      # Validate
      if (nchar(guide_seq_rc_chr) != 20 || nchar(pam_seq_rc_chr) != 3 ||
          substr(pam_seq_rc_chr, 2, 3) != "GG" || nchar(full_seq_rc_chr) != 23) {
        next
      }

      # Convert RC coordinates to original coordinates
      pam_end_orig_int <- search_start_pos_int + search_seq_length_bp_int - pam_start_rc_int
      pam_start_orig_int <- search_start_pos_int + search_seq_length_bp_int - pam_end_rc_int
      guide_end_orig_int <- search_start_pos_int + search_seq_length_bp_int - guide_start_rc_int
      guide_start_orig_int <- search_start_pos_int + search_seq_length_bp_int - guide_end_rc_int

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

  # Combine guides
  gene_guides_lst <- c(guides_fwd_lst, guides_rev_lst)

  cat("  Guides extracted:", length(gene_guides_lst), "\n")

  # ==========================================================================
  # CDS overlap and synonymous mutation analysis
  # ==========================================================================

  if (length(gene_guides_lst) > 0) {
    # Extract CDS sequence once
    cds_seq_dna <- Biostrings::subseq(x = seq_dna,
                                       start = CDS_START_EXPECTED_pos_int,
                                       end = cds_end_pos_int)
    cds_seq_chr <- as.character(cds_seq_dna)

    for (guide_idx in seq_along(gene_guides_lst)) {
      guide_info_lst <- gene_guides_lst[[guide_idx]]

      # Check CDS overlap
      guide_region_start_int <- guide_info_lst$guide_start_pos
      guide_region_end_int <- guide_info_lst$pam_end_pos

      overlaps_cds_lgl <- (guide_region_start_int <= cds_end_pos_int &&
                           guide_region_end_int >= CDS_START_EXPECTED_pos_int)

      guide_info_lst$overlaps_cds <- overlaps_cds_lgl

      if (overlaps_cds_lgl) {
        guide_info_lst$cds_start <- CDS_START_EXPECTED_pos_int
        guide_info_lst$cds_end <- cds_end_pos_int

        # Determine which PAM positions to check
        if (guide_info_lst$guide_strand == "+") {
          pam_check_positions_int <- 2:3
          pam_for_analysis_chr <- guide_info_lst$pam_sequence
        } else {
          pam_check_positions_int <- 1:2
          pam_for_analysis_chr <- as.character(
            Biostrings::reverseComplement(
              Biostrings::DNAString(guide_info_lst$pam_sequence)
            )
          )
        }

        # Check PAM mutability
        pam_positions_int <- guide_info_lst$pam_start_pos:guide_info_lst$pam_end_pos
        pam_can_mutate_lgl <- FALSE

        for (pam_idx in pam_check_positions_int) {
          pam_abs_pos_int <- pam_positions_int[pam_idx]
          pos_in_cds_int <- pam_abs_pos_int - CDS_START_EXPECTED_pos_int + 1

          if (pos_in_cds_int >= 1 && pos_in_cds_int <= nchar(cds_seq_chr)) {
            codon_number_int <- ceiling(pos_in_cds_int / 3)
            position_in_codon_int <- ((pos_in_cds_int - 1) %% 3) + 1

            codon_start_int <- (codon_number_int - 1) * 3 + 1
            codon_end_int <- codon_start_int + 2

            if (codon_end_int <= nchar(cds_seq_chr)) {
              codon_chr <- substr(cds_seq_chr, codon_start_int, codon_end_int)
              current_aa_chr <- GENETIC_CODE_nls[[codon_chr]]

              if (!is.null(current_aa_chr) && !is.na(current_aa_chr)) {
                current_nt_chr <- substr(codon_chr, position_in_codon_int, position_in_codon_int)

                for (new_nt_chr in c("A", "T", "G", "C")) {
                  if (new_nt_chr != current_nt_chr) {
                    mutated_codon_chr <- codon_chr
                    substr(mutated_codon_chr, position_in_codon_int, position_in_codon_int) <- new_nt_chr
                    mutated_aa_chr <- GENETIC_CODE_nls[[mutated_codon_chr]]

                    if (!is.null(mutated_aa_chr) && !is.na(mutated_aa_chr) &&
                        mutated_aa_chr == current_aa_chr) {
                      pam_can_mutate_lgl <- TRUE
                      break
                    }
                  }
                }
              }
            }
          }
          if (pam_can_mutate_lgl) break
        }

        guide_info_lst$pam_gg_can_mutate <- pam_can_mutate_lgl

        # Count wobble positions
        guide_positions_int <- guide_info_lst$guide_start_pos:guide_info_lst$guide_end_pos
        guide_in_cds_int <- guide_positions_int[
          guide_positions_int >= CDS_START_EXPECTED_pos_int &
          guide_positions_int <= cds_end_pos_int
        ]

        if (length(guide_in_cds_int) > 0) {
          pos_in_cds_int <- guide_in_cds_int - CDS_START_EXPECTED_pos_int + 1
          position_in_codon_int <- ((pos_in_cds_int - 1) %% 3) + 1
          wobble_count_int <- sum(position_in_codon_int == 3)
        } else {
          wobble_count_int <- 0
        }

        guide_info_lst$guide_wobble_sites <- wobble_count_int

        # Determine difficulty
        if (guide_info_lst$pam_gg_can_mutate) {
          guide_info_lst$silencing_difficulty <- "easy"
          guide_info_lst$silencing_strategy <- "PAM_mutation"
        } else if (wobble_count_int >= 3) {
          guide_info_lst$silencing_difficulty <- "moderate"
          guide_info_lst$silencing_strategy <- "guide_mutations"
        } else if (wobble_count_int >= 1) {
          guide_info_lst$silencing_difficulty <- "difficult"
          guide_info_lst$silencing_strategy <- "limited_guide_mutations"
        } else {
          guide_info_lst$silencing_difficulty <- "very_difficult"
          guide_info_lst$silencing_strategy <- "non_synonymous_required"
        }

      } else {
        # Non-coding
        guide_info_lst$cds_start <- NA
        guide_info_lst$cds_end <- NA
        guide_info_lst$pam_gg_can_mutate <- TRUE
        guide_info_lst$guide_wobble_sites <- NA
        guide_info_lst$silencing_difficulty <- "easy"
        guide_info_lst$silencing_strategy <- "non_coding_region"
      }

      gene_guides_lst[[guide_idx]] <- guide_info_lst
    }

    # Store results for this gene
    all_genes_results_lst[[gene_name_chr]] <- gene_guides_lst

    cat("  CDS overlap:", sum(sapply(gene_guides_lst, function(x) x$overlaps_cds)), "/",
        length(gene_guides_lst), "\n")
  }
}

cat("\n=== ALL FILES PROCESSED ===\n")
cat("Total genes processed:", length(all_genes_results_lst), "\n")
# ============================================================================
# STEP 15: Convert results to data frame and write output
# ============================================================================

cat("\n=== CREATING OUTPUT FILE ===\n")

# Flatten all results into single list
all_guides_flat_lst <- unlist(all_genes_results_lst, recursive = FALSE)

cat("Total guides across all genes:", length(all_guides_flat_lst), "\n")

# Convert to data frame
results_df <- data.frame(
  gene_name = sapply(all_guides_flat_lst, function(x) x$gene_name),
  guide_strand = sapply(all_guides_flat_lst, function(x) x$guide_strand),
  guide_sequence = sapply(all_guides_flat_lst, function(x) x$guide_sequence),
  pam_sequence = sapply(all_guides_flat_lst, function(x) x$pam_sequence),
  full_sequence = sapply(all_guides_flat_lst, function(x) x$full_sequence),
  guide_start_pos = sapply(all_guides_flat_lst, function(x) x$guide_start_pos),
  guide_end_pos = sapply(all_guides_flat_lst, function(x) x$guide_end_pos),
  pam_start_pos = sapply(all_guides_flat_lst, function(x) x$pam_start_pos),
  pam_end_pos = sapply(all_guides_flat_lst, function(x) x$pam_end_pos),
  overlaps_cds = sapply(all_guides_flat_lst, function(x) x$overlaps_cds),
  cds_start = sapply(all_guides_flat_lst, function(x) 
    ifelse(is.null(x$cds_start) || is.na(x$cds_start), NA, x$cds_start)),
  cds_end = sapply(all_guides_flat_lst, function(x) 
    ifelse(is.null(x$cds_end) || is.na(x$cds_end), NA, x$cds_end)),
  pam_gg_can_mutate = sapply(all_guides_flat_lst, function(x) x$pam_gg_can_mutate),
  guide_wobble_sites = sapply(all_guides_flat_lst, function(x) 
    ifelse(is.null(x$guide_wobble_sites) || is.na(x$guide_wobble_sites), NA, x$guide_wobble_sites)),
  silencing_difficulty = sapply(all_guides_flat_lst, function(x) x$silencing_difficulty),
  silencing_strategy = sapply(all_guides_flat_lst, function(x) x$silencing_strategy),
  stringsAsFactors = FALSE
)

# Show summary
cat("\nData frame created with", nrow(results_df), "rows and", ncol(results_df), "columns\n")

cat("\nFinal summary statistics:\n")
cat("Genes analyzed:", length(unique(results_df$gene_name)), "\n")
cat("Forward strand guides:", sum(results_df$guide_strand == "+"), "\n")
cat("Reverse strand guides:", sum(results_df$guide_strand == "-"), "\n")
cat("Guides overlapping CDS:", sum(results_df$overlaps_cds), "\n")
cat("Guides in non-coding regions:", sum(!results_df$overlaps_cds), "\n")

cat("\nSilencing difficulty distribution:\n")
print(table(results_df$silencing_difficulty))

# Write to TSV file
cat("\nWriting output to:", OUTPUT_FILE_PATH_chr, "\n")

utils::write.table(x = results_df,
                   file = OUTPUT_FILE_PATH_chr,
                   sep = "\t",
                   row.names = FALSE,
                   quote = FALSE,
                   na = "NA")

cat("\n=== OUTPUT COMPLETE ===\n")
cat("File written:", OUTPUT_FILE_PATH_chr, "\n")
cat("Total guides:", nrow(results_df), "\n")

# Show first few rows as preview
cat("\nPreview of output (first 5 rows):\n")
print(head(results_df, 5))
