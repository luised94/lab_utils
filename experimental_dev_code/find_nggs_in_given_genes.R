################################################################################
# Find NGG PAM Sites in Target Genes for S. cerevisiae
# Author: Luis | Date: 2025-11-11 | Version: 1.0.0
################################################################################
# PURPOSE:
#   Find NGG in target genes given a reference genome and gff annotation file.
# USAGE:
#   Source from repl.
# DEPENDENCIES:
#   See libraries, requires genome and gff file.
# OUTPUTS:
#   See all_pam_sites to look at results.
################################################################################

standard_name_field <- "gene" # or Alias
systematic_name_field <- "ID"
TYPE_OF_FEATURE <- "gene"
# Define the search pattern.
SEARCH_PATTERN <- "NGG" # NGG for PAM sites in Cas9.
# --- 1. Load Required Libraries ---
# These packages are the standard for genomic analysis in R.
library(rtracklayer)
library(Biostrings)
library(GenomicRanges)

# --- 2. USER INPUT: Define Paths and Gene List ---
# Define the file paths for your reference genome and GFF annotation file.
GFF_FILEPATH <- "/home/luised94/data/feature_files/240830_saccharomyces_cerevisiae.gff"
GENOME_FILEPATH <- "/home/luised94/data/REFGENS/SaccharomycescerevisiaeS288C/SaccharomycescerevisiaeS288C_refgenome.fna"

# Define the list of genes you want to analyze.
# Colloquial gene names
TARGET_GENES <- c(
  "Toa1", "Sua7",
  "Toa2", "Spt15"
)

#possible_fields <- c("Name", "gene_name", "gene", "Alias")


# --- 3. Load Genome and Annotation Files ---
# Check if the specified files exist before running the script.
if (!file.exists(GFF_FILEPATH) || !file.exists(GENOME_FILEPATH)) {
  stop("Error: GFF or Genome file not found. Please verify your file paths.")
}

cat("Loading GFF annotation file...\n")
gff_data <- rtracklayer::import(GFF_FILEPATH)
cat("Loading reference genome FASTA file...\n")
genome_seq <- Biostrings::readDNAStringSet(GENOME_FILEPATH)


# Clean up the names of the genome sequences (chromosomes).
names(genome_seq) <- sapply(strsplit(names(genome_seq), " "), `[`, 1)

# --- 4. Main Loop to Process Each Gene ---

# This list will hold the final results for all processed genes.
all_pam_sites <- list()

cat("\nStarting to process genes...\n")
genes_data <- gff_data[
  gff_data$type == "gene" &
  !is.na(gff_data$type == "gene")
]

for (gene_name in TARGET_GENES) {

  cat(paste0("  -> Processing gene: ", gene_name, "\n"))

  # --- A. Find the gene's coordinates in the GFF file ---
  # Filter the GFF data for an entry with type 'gene' that matches the current gene name.
  # Note: The attribute containing the gene name in a GFF can vary (e.g., 'Name', 'ID', 'gene_id').

  col_data <- mcols(genes_data)[[standard_name_field]]
  matches <- grepl(gene_name, col_data, ignore.case = TRUE)
  # Handle CharacterList columns (which return LogicalList)
  if (is(matches, "LogicalList")) {
    matches <- any(matches)
  }

  gene_coords <- genes_data[matches]

  if (length(gene_coords) == 0) {
    warning(paste("Could not find coordinates for gene:", gene_name, ". Skipping."))
    next # Skip to the next gene in the list
  }

  if (length(gene_coords) > 1) {
    warning(paste("More than one match found for gene name:", gene_name, ". Taking first match."))
  }

  # Use the first match if multiple are found
  gene_coords <- gene_coords[1, ]

  # Get the chromosome name from the gene's coordinates.
  seq_name <- as.character(GenomicRanges::seqnames(gene_coords))

  # Verify that this chromosome exists in our loaded genome file.
  if (!seq_name %in% names(genome_seq)) {
      warning(paste("Chromosome '", seq_name, "' for gene '", gene_name, "' not found in genome file. Skipping."))
      next
  }

  # --- B. Extract the gene's DNA sequence ---


  # --- C. Find all "NGG" PAM sites ---
 gene_sequence <- subseq(genome_seq[[seq_name]],
                         start(gene_coords),
                         end(gene_coords))

  gene_sequence <- Biostrings::getSeq(genome_seq, gene_coords)
  pam_matches <- Biostrings::matchPattern(SEARCH_PATTERN, gene_sequence, fixed = FALSE)

  # --- D. Store the results ---

  if (length(pam_matches[[1]]) > 0) {
    # If matches were found, create a clean data frame with the results.
    gene_results_df <- data.frame(
      Gene = gene_name,
      PAM_Start_in_Gene = Biostrings::start(pam_matches[[1]]),
      PAM_End_in_Gene = Biostrings::end(pam_matches[[1]]),
      PAM_Sequence = as.character(
        Biostrings::DNAStringSet(
          gene_sequence,
          start =  Biostrings::start(pam_matches[[1]]),
          width = 3
        )
      )
      #PAM_Sequence = as.character(
      #    Biostrings::subseq(
      #      gene_sequence,
      #      Biostrings::start(pam_matches[[1]]),
      #      Biostrings::end(pam_matches[[1]]
      #)))
    )
    # Add the results to our main list, named by the gene.
    all_pam_sites[[gene_name]] <- gene_results_df
  } else {
    # If no matches found, store a NULL or empty data frame.
    all_pam_sites[[gene_name]] <- data.frame()
  }

} # End of the for-loop

# --- 5. Display the Final Results ---
cat("\n--- PAM Site Search Complete ---\n\n")

if (length(all_pam_sites) > 0) {
  # The `print()` function will display the list of data frames in a readable format.
  print(all_pam_sites)
} else {
  cat("No results were generated. Please check your gene list and file paths.\n")
}
