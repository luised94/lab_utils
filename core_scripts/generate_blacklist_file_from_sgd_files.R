# @USAGE source("core_scripts/generate_blacklist_file_from_sgd_files.R")
#---------------------------------------------------------------------
# CONFIGURATION AND CONSTANTS
#---------------------------------------------------------------------
# Current date for output filename
current_date <- format(Sys.time(), "%Y%m%d")

# File directory
file_directory <- normalizePath("~/data/feature_files")

# Output file path
output_file_path <- file.path(file_directory, paste0(current_date, "_saccharomyces_cerevisiae_s288c_blacklist.bed"))

blacklist_patterns <- c(
  # Repetitive elements
  "long_terminal_repeat",
  "ltr_retrotransposon",
  "telomeric_repeat",
  "x_element_combinatorial_repeat",  # Note: Order matters - put this before x_element to avoid partial matches
  "x_element",
  "y_prime_element",
  "transposable_element_gene",
  # Special genomic features
  "telomere",
  "centromere",
  "mating_type_region",
  "silent_mating_type_cassette_array",
  # Specialized regions
  "w_region",
  "x_region",
  "y_region",
  "z1_region",
  "z2_region",
  # RNA genes and related regions - highly important for sponge sequences
  "rrna_gene",
  "trna_gene",
  "snrna_gene",
  "snorna_gene",
  "internal_transcribed_spacer_region",
  "external_transcribed_spacer_region",
  "non_transcribed_region",
  # Other potentially problematic features
  "telomerase_rna_gene",
  "centromere_dna_element_i",
  "centromere_dna_element_ii",
  "centromere_dna_element_iii"
)

#---------------------------------------------------------------------
# PREPROCESSING AND SETUP
#---------------------------------------------------------------------

# Get all .bed files in the directory
all_files <- list.files(path = file_directory, pattern = "\\.bed$", full.names = FALSE)

# Filter files to include only those matching our blacklist patterns
# Each file should match the pattern: [timestamp]_[pattern]_sgd.bed
blacklist_files <- character(0)

# Track which patterns were found
patterns_found <- logical(length(blacklist_patterns))
names(patterns_found) <- blacklist_patterns

# Build the list of files to process
for (pattern_index in 1:length(blacklist_patterns)) {
  pattern <- blacklist_patterns[pattern_index]

  # Create regex that matches timestamp_[pattern]_sgd.bed
  regex_pattern <- paste0("^\\d+_", pattern, "_sgd\\.bed$")

  # Find matching files
  matching_files <- all_files[grep(regex_pattern, all_files)]

  if (length(matching_files) > 0) {
    blacklist_files <- c(blacklist_files, matching_files)
    patterns_found[pattern_index] <- TRUE

    # Report matches
    cat("Found", length(matching_files), "files matching pattern:", pattern, "\n")
  } else {
    cat("WARNING: No files found matching pattern:", pattern, "\n")
  }
} # end for loop for finding files

# Report patterns not found
missing_patterns <- blacklist_patterns[!patterns_found]
if (length(missing_patterns) > 0) {
  cat("WARNING: No files found for the following patterns:\n")
  cat(paste(" -", missing_patterns), sep = "\n")
}

# Check if we found any files to process
if (length(blacklist_files) == 0) {
  stop("No matching blacklist files found in directory: ", file_directory)
}

# Convert relative paths to absolute paths using the directory
blacklist_files_full_paths <- file.path(file_directory, blacklist_files)

# Check if all files exist
missing_files <- blacklist_files_full_paths[!file.exists(blacklist_files_full_paths)]
if (length(missing_files) > 0) {
  stop("The following files are missing: ", paste(missing_files, collapse = ", "))
}

# Initialize storage for combined bed records
combined_bed_data <- data.frame(
  chromosome = character(),
  start = numeric(),
  end = numeric(),
  name = character(),
  score = numeric(),
  strand = character(),
  stringsAsFactors = FALSE
)

# Initialize counters for reporting
total_records_read <- 0
files_processed <- 0

cat("Starting to process", length(blacklist_files), "blacklist files\n")

#---------------------------------------------------------------------
# MAIN PROCESSING LOGIC
#---------------------------------------------------------------------

# Read and combine each blacklist file
for (bed_file_index in 1:length(blacklist_files)) {
  bed_file_path <- blacklist_files_full_paths[bed_file_index]
  bed_file_name <- blacklist_files[bed_file_index]

  # Report progress
  cat("Processing file:", bed_file_name, "\n")

  # Try to read the file with error handling
  tryCatch({
    # Read BED file - handling variable number of columns
    # Standard BED format has at least 3 columns: chromosome, start, end
    # May also have name, score, strand, etc.
    bed_content <- read.table(bed_file_path, header = FALSE, sep = "\t", 
                              stringsAsFactors = FALSE, comment.char = "#",
                              fill = TRUE, quote = "")

  # Filter out rows where start or end positions aren't numeric
  #valid_rows <- !is.na(as.numeric(bed_content[,2])) & !is.na(as.numeric(bed_content[,3]))
  #bed_content <- bed_content[valid_rows, ]
    # Ensure we have at least the required 3 columns
    if (ncol(bed_content) < 3) {
      warning("File ", bed_file_name, " has fewer than 3 columns. Skipping.")
      next
    }

    # # Standardize column names for the first 6 BED columns
    # valid_columns <- min(ncol(bed_content), 6)
    # names(bed_content)[1:valid_columns] <- c("chromosome", "start", "end", 
    #                                         "name", "score", "strand")[1:valid_columns]
    # 
    # # If name column doesn't exist, create it using the file basename
    # if (valid_columns < 4) {
    #   bed_content$name <- rep(sub("\\.bed$", "", basename(bed_file_path)), nrow(bed_content))
    #   valid_columns <- 4  # Now we have 4 columns
    # }
    # 
    # # If score column doesn't exist, set default as 0
    # if (valid_columns < 5) {
    #   bed_content$score <- 0
    #   valid_columns <- 5  # Now we have 5 columns
    # }
    # 
    # # If strand column doesn't exist, set default as "."
    # if (valid_columns < 6) {
    #   bed_content$strand <- "."
    #   valid_columns <- 6  # Now we have 6 columns
    # }

    # Assume the input files are standardized as mentioned
    # Just use the first 3 columns which are the essential BED format
    standardized_data <- data.frame(
      chromosome = bed_content[,1],
      start = as.numeric(bed_content[,2]),
      end = as.numeric(bed_content[,3]),
      name = if(ncol(bed_content) >= 4) bed_content[,4] else rep(sub("\\.bed$", "", basename(bed_file_name)), nrow(bed_content)),
      score = if(ncol(bed_content) >= 5) bed_content[,5] else 0,
      strand = if(ncol(bed_content) >= 6) bed_content[,6] else ".",
      feature_type = if(ncol(bed_content) >= 7) bed_content[,7] else ".",
      gene_name = if(ncol(bed_content) >= 8) bed_content[,8] else ".",
      stringsAsFactors = FALSE
    )

    # Update records count
    records_in_file <- nrow(standardized_data)
    total_records_read <- total_records_read + records_in_file

    # Append to the combined dataset
    combined_bed_data <- rbind(combined_bed_data, standardized_data)

    # Update processed files counter
    files_processed <- files_processed + 1

    # Report for this file
    cat("  Added", records_in_file, "records from", bed_file_name, "\n")

  }, error = function(e) {
    warning("Error processing file ", bed_file_name, ": ", e$message)
  })
}

#---------------------------------------------------------------------
# DATA VALIDATION AND ASSERTIONS
#---------------------------------------------------------------------
# Verify we processed some files
if (files_processed == 0) {
  stop("No files were successfully processed. Check file paths and formats.")
}

# Verify we have data
if (nrow(combined_bed_data) == 0) {
  stop("No records were read from the input files.")
}

# Check for any invalid BED format entries (end should be >= start)
invalid_ranges <- which(combined_bed_data$end < combined_bed_data$start)
if (length(invalid_ranges) > 0) {
  warning("Found ", length(invalid_ranges), " invalid ranges (end < start). Fixing...")
  # Swap start and end for invalid ranges
  for (idx in invalid_ranges) {
    temp <- combined_bed_data$start[idx]
    combined_bed_data$start[idx] <- combined_bed_data$end[idx]
    combined_bed_data$end[idx] <- temp
  }
}

#---------------------------------------------------------------------
# SORT BED FILE
#---------------------------------------------------------------------
# Sort BED entries by chromosome and start position
# First convert chromosome names to factors to ensure proper sorting
# when chromosomes are named like "chr1", "chr2", etc.
combined_bed_data$sortkey <- gsub("chr", "", combined_bed_data$chromosome)
combined_bed_data$sortkey <- as.numeric(ifelse(grepl("^\\d+$", combined_bed_data$sortkey), 
                                              combined_bed_data$sortkey, 999))
combined_bed_data <- combined_bed_data[order(combined_bed_data$sortkey, combined_bed_data$start), ]
combined_bed_data$sortkey <- NULL  # Remove the temporary sort key
#combined_bed_data <- combined_bed_data[!is.na(combined_bed_data$start) & !is.na(combined_bed_data$end), ]

#---------------------------------------------------------------------
# REMOVE REDUNDANT SECTIONS
#---------------------------------------------------------------------
library(GenomicRanges)
library(rtracklayer)

# Step 1: Report the size of the original dataset
cat("Processing", nrow(combined_bed_data), "total regions from all files\n")
message("Removing redundant regions...")
cat("Converting to GRanges object...\n")
gr <- GenomicRanges::GRanges(
  seqnames = combined_bed_data$chromosome,
  ranges = IRanges::IRanges(
    start = combined_bed_data$start + 1,  # BED is 0-based, GRanges is 1-based
    end = combined_bed_data$end
  ),
  strand = combined_bed_data$strand,
  # Include all other columns as metadata
  name = combined_bed_data$name,
  score = combined_bed_data$score,
  feature_type = combined_bed_data$feature_type,
  gene_name = combined_bed_data$gene_name
)

cat("Removing redundancy by merging overlapping regions...\n")
non_redundant_gr <- GenomicRanges::reduce(gr, with.revmap=TRUE)

# Step 5: Report the reduction
cat("Original regions:", length(gr), "\n")
cat("Non-redundant regions:", length(non_redundant_gr), "\n")
cat("Reduction:", round((1 - length(non_redundant_gr)/length(gr)) * 100, 2), "%\n")

# Step: 4: Add metadata back to the reduced regions
new_mcols <- S4Vectors::DataFrame(
  name = character(length(non_redundant_gr)),
  score = numeric(length(non_redundant_gr)),
  feature_type = character(length(non_redundant_gr)),
  gene_name = character(length(non_redundant_gr))
)

# For each reduced range, aggregate metadata from the original ranges that map to it
for (i in seq_along(non_redundant_gr)) {
  # Get indices of original ranges that were merged into this reduced range
  idx <- non_redundant_gr$revmap[[i]]
  for (mcols in names(new_mcols)) {
    if (mcols == "score"){
      new_mcols[i, mcols] <- max(GenomicRanges::mcols(gr)[idx, mcols])
    } else {
      x <- na.omit(GenomicRanges::mcols(gr)[idx, mcols])
      if (length(x) == 1 || all(x == x[1])) {
        # If all values are the same, just return one
        column_metadata <- x[1]
      } else {
        # Get unique values and collapse with commas
        column_metadata <- paste(unique(x), collapse=",")
      }
      new_mcols[i, mcols] <- column_metadata
    }
  }
}

# Assign the new metadata to the reduced GRanges
GenomicRanges::mcols(non_redundant_gr) <- new_mcols

# Remove the revmap column which we no longer need
non_redundant_gr$revmap <- NULL

#---------------------------------------------------------------------
# WRITE OUTPUT FILES
#---------------------------------------------------------------------
# When exporting to BED, ensure the column order matches your expectations
# Export the non-redundant regions as BED file
cat("Exporting blacklist BED file...\n")
#rtracklayer::export(non_redundant_gr, output_file_path, format = "BED")
# First create the bed_data data frame
non_redundant_bed <- data.frame(
  chromosome = as.character(GenomicRanges::seqnames(non_redundant_gr)),
  start = GenomicRanges::start(non_redundant_gr) - 1,  # Convert back to 0-based
  end = GenomicRanges::end(non_redundant_gr),
  name = GenomicRanges::mcols(non_redundant_gr)$name,
  score = GenomicRanges::mcols(non_redundant_gr)$score,
  strand = as.character(GenomicRanges::strand(non_redundant_gr)),
  feature_type = GenomicRanges::mcols(non_redundant_gr)$feature_type,
  gene_name = GenomicRanges::mcols(non_redundant_gr)$gene_name,
  stringsAsFactors = FALSE
)

# Write the data frame to a BED file
write.table(
  non_redundant_bed,
  file = output_file_path,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

cat("Exported BED file with all columns preserved\n")

# For documentation purposes, save the column names to a separate file
#write.table(
#  data.frame(column_name = names(non_redundant_bed)),
#  file = "yeast_blacklist_sponge.columns.txt",
#  quote = FALSE,
#  sep = "\t",
#  row.names = FALSE,
#  col.names = FALSE
#)
#
#cat("Column names saved to yeast_blacklist_sponge.columns.txt\n")

# Write the combined data to output file
#cat("Writing", nrow(combined_bed_data), "records to output file:", output_file_path, "\n")
#write.table(combined_bed_data, file = output_file_path, 
#            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Final verification
if (file.exists(output_file_path)) {
  file_info <- file.info(output_file_path)
  cat("Successfully created blacklist file:", output_file_path, 
      "\nFile size:", file_info$size, "bytes\n",
      "Total records written:", nrow(combined_bed_data), "\n",
      "Files processed:", files_processed, "of", length(blacklist_files), "\n")
} else {
  stop("Failed to create output file:", output_file_path)
}

cat("Blacklist generation complete.\n")
cat("Run bedtools merge -i my_blacklist_original.bed > my_blacklist_merged.bed to ensure bed file is merged.")
