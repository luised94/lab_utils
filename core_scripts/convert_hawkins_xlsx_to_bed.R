###############################################################################
# Turn hawkins data set into bed file
################################################################################
# PURPOSE: Output hawkins dataset as bed file to have as reference
# Conclusion: {{fill}}
# USAGE: source("reference_code/pipeline_completion/plot_bigwig_files.R")
# DEPENDENCIES: 
# packages = GenomicRanges, rtracklayer
# data = hawkins 2013 excel sheet download using curl
# OUTPUT: Bed files for different gene feature sets and the sgd file per se
# AUTHOR: LEMR
# DATE: 2025-04-07
################################################################################
# Library start ---------------------------------------------------
library(tidyverse)
library(dplyr)
# [2] Main Processing Script ---------------------------------------------------
input_file_path <- file.path(Sys.getenv("HOME"), "data", "feature_files", "240830_hawkins-origins-timing.xlsx")

# Initial checks
stopifnot(
  "Input file does not exist" = file.exists(input_file_path),
  "Input file is not readable" = file.access(input_file_path, 4) == 0
)

cat("Reading input file with quote handling...\n")
tryCatch({
  raw_data <- readxl::read_excel(input_file_path)
}, error = function(e) {
  cat("FILE READING FAILED:\n")
  cat("First 5 lines from readLines():\n")
  print(head(readLines(input_file_path), 5))
  stop(e)
})

hawkins_df <- as.data.frame(raw_data) %>%
    dplyr::select(1:7) %>%
    filter(!is.na(Chromosome))

cat(sprintf("Detected %d columns\n", ncol(raw_data)))
cat(sprintf("Successfully read %d rows with %d columns\n", 
            nrow(hawkins_df), ncol(hawkins_df)))
stopifnot(
  "Unexpected column count (expected 7)" = ncol(hawkins_df) == 7
)

required_bed_fields <- c("chrom", "start", "end", "name", "score", "strand", "feature_type")
hawkins_df %>% mutate()
#Construct a GRange object from timing dataframe
#GRanges(seqnames = timing$Chromosome, ranges = IRanges(start = timing$Position-100, end = timing$Position+100), strand = "*", timing$T1/2)
#
#
## Column processing with validation
#raw_data$chromosome_number <- as.integer(raw_data[[9]])
#raw_data$start <- as.integer(raw_data[[10]])
#raw_data$end <- as.integer(raw_data[[11]])
#raw_data$strand_code <- raw_data[[12]]
#feature_types <- raw_data[[2]]
#
## Validation checks
#stopifnot(
#  "Missing chromosome numbers" = !any(is.na(raw_data$chromosome_number)),
#  "Invalid start coordinates" = !any(is.na(raw_data$start)),
#  "Invalid end coordinates" = !any(is.na(raw_data$end)),
#  "Invalid strand codes" = all(raw_data$strand_code %in% c("W", "C"))
#)
#
## Chromosome conversion to roman numerals
#valid_chromosomes <- c(1:17)
#invalid_chrom <- unique(raw_data$chromosome_number[!raw_data$chromosome_number %in% valid_chromosomes])
#stopifnot(
#  "Invalid chromosome numbers detected" = length(invalid_chrom) == 0
#)
#cat("Converting chromosome numbers to roman numerals...\n")
#raw_data$chrom <- paste0(
#  "chr", 
#  as.roman(raw_data$chromosome_number)
#)
#
#cat("Processing genomic coordinates...\n")
#raw_data$bed_start <- pmin(raw_data$start, raw_data$end) - 1L
#raw_data$bed_end <- pmax(raw_data$start, raw_data$end)
#raw_data$bed_strand <- ifelse(raw_data$strand_code == "W", "+", "-")
#
#cat("Resolving feature names...\n")
#
## Create ORF name lookup based on shared identifiers (V1)
#orf_name_lookup <- raw_data[
#  raw_data$V2 == "ORF" & !is.na(raw_data$V4),
#  c("V1", "V4")
#]
#colnames(orf_name_lookup) <- c("source_id", "canonical_name")
#
#merged_data <- merge(
#  raw_data,
#  orf_name_lookup,
#  by.x = "V1",
#  by.y = "source_id",
#  all.x = TRUE
#)
#
#merged_data$resolved_name <- dplyr::coalesce(
#  merged_data$V4,
#  merged_data$canonical_name,
#  # S000350094_CDS format
#  paste0(merged_data$V1, "_", merged_data$V2)
#)
#
#bed_data <- data.frame(
#    chrom = merged_data$chrom,
#    start = merged_data$bed_start,
#    end = merged_data$bed_end,
#    name = merged_data$resolved_name,
#    score = 0,
#    strand = merged_data$bed_strand,
#    feature_type = merged_data$V2,
#    gene_name = merged_data$V5,
#    stringsAsFactors = FALSE
#)
#
#stopifnot(
#  "Invalid BED coordinates (start >= end)" = all(bed_data$start < bed_data$end),
#  "Missing chromosome names" = !any(is.na(bed_data$chrom)),
#  "Empty feature names" = !any(is.na(bed_data$name) | bed_data$name == "")
#)
#
#cat("Checking for duplicate features...\n")
#duplicate_keys <- paste(bed_data$chrom, bed_data$start, bed_data$end, bed_data$name)
#if(any(duplicated(duplicate_keys))) {
#  dup_count <- sum(duplicated(duplicate_keys))
#  cat(sprintf("WARNING: Found %d duplicate features\n", dup_count))
#  print(head(bed_data[duplicated(duplicate_keys), ]))
#  cat("Consider adding feature type suffixes\n")
#  
#  # Add feature type disambiguation
#  bed_data$name <- paste0(bed_data$name, "_", merged_data$V2)
#}
#
#cat("\nWriting output files:\n")
#unique_types <- unique(feature_types)
#cat(sprintf("Found %d feature types: %s\n",
#            length(unique_types), paste(unique_types, collapse = ", ")))
#
#output_directory <- normalizePath(
#  file.path("~", "data", "feature_files"),
#  mustWork = FALSE
#)
#dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)
#
## Verify directory creation
#stopifnot(
#  "Output directory creation failed" = dir.exists(output_directory),
#  "Output directory not writable" = file.access(output_directory, 2) == 0
#)
#
## Track base names and their counts
#name_registry <- list()
#
## Generate manifest with proper numbering
#manifest <- data.frame()
#current_date <- format(Sys.Date(), "%Y%m%d")
#
#for (feature_type in unique(feature_types)) {
#
#    # Sanitize file name
#  # Convert to lowercase and replace special characters
#  clean_name <- tolower(feature_type)
#  # Replace all non-alphanumeric characters with underscores
#  clean_name <- gsub("[^a-z0-9_]+", "_", clean_name)
#  # Collapse multiple underscores
#  clean_name <- gsub("_+", "_", clean_name)
#  # Trim leading/trailing underscores
#  base_name <- gsub("^_|_$", "", clean_name)
#  # Update registry to track file names
#  if (is.null(name_registry[[base_name]])) {
#    name_registry[[base_name]] <- 1
#    file_name <- paste0(current_date, "_", base_name, "_sgd.bed")
#  } else {
#    count <- name_registry[[base_name]]
#    name_registry[[base_name]] <- count + 1
#    file_name <- paste0(current_date, "_", base_name, "_", count, "_sgd.bed")
#  }
#  count <- sum(feature_types == feature_type)
#
#  manifest <- rbind(manifest, data.frame(
#    Original_Feature = feature_type,
#    Sanitized_Name = file_name,
#    Entries = count
#  ))
#}
#
## Show final manifest
#cat("\nFinal File Manifest:\n")
#print(manifest)
#
## Interactive verification
#response <- tolower(readline(prompt = "Proceed with writing files? (y/n): "))
#if (!response %in% c("y", "yes")) {
#  cat("Operation cancelled\n")
#    stop("Inspect files or resource the script.")
#  #quit(save = "no", status = 0)
#}
#
## Actual file writing
#cat("\nWriting files...\n")
#for (i in 1:nrow(manifest)) {
#  current <- manifest[i,]
#  type_mask <- feature_types == current$Original_Feature
#  current$FilePath <- file.path(output_directory, current$Sanitized_Name)
#  write.table(
#    x = bed_data[type_mask, ],
#    file = current$FilePath,
#    sep = "\t",
#    quote = FALSE,
#    row.names = FALSE,
#    col.names = TRUE
#  )
#
#  cat(sprintf(" - %-25s: %6d entries -> %s\n",
#              current$Original_Feature,
#              current$Entries,
#              current$FilePath))
#}
#
#cat(sprintf("\nSuccessfully wrote %d files to:\n%s\n",
#            nrow(manifest), output_directory))
#
#cat("\n=== COMBINED FILE GENERATION ===\n")
#combined_path <- file.path(output_directory,
#    paste(current_date, "_", "sgd_features_combined.bed", sep = "")
#)
#
## Interactive confirmation
#cat(sprintf(
#  "Propose creating combined file with %d features:\n%s\n",
#  nrow(bed_data),
#  combined_path
#))
#
#create_combined <- tolower(readline("Create combined BED file? (y/n): "))
#if(create_combined %in% c("y", "yes")) {
#  ## Write with track line for genome browsers
#  #track_line <- sprintf(
#  #  'track name="SGD_Features" description="Combined_SGD_Annotation" visibility=2 itemRgb="On"'
#  #)
#
#  # Prepare final structure with Gviz-friendly columns
#  #combined_output <- bed_data[, c("chrom", "start", "end", "name", "score", "strand", "feature_type")]
#  combined_output <- bed_data
#  #combined_output$thickStart <- combined_output$start
#  #combined_output$thickEnd <- combined_output$end
#  #combined_output$itemRgb <- ifelse(
#  #  bed_data$feature_type == "ORF",
#  #  "255,0,0",  # Red for ORFs
#  #  "0,0,255"    # Blue for others
#  #)
#
#  # Write to file
#  #writeLines(track_line, combined_path)
#  write.table(
#    combined_output,
#    file = combined_path,
#    append = TRUE,
#    sep = "\t",
#    quote = FALSE,
#    row.names = FALSE,
#    col.names = TRUE
#  )
#  
#  cat(sprintf(
#    "\nSuccessfully wrote combined BED file with:\n- %d total features\n- Color-coded types\n- Genome browser compatibility\n",
#    nrow(combined_output)
#  ))
#  
#  # Create matching Gviz style file
#  style_path <- file.path(output_directory, "gviz_style.txt")
#  writeLines(
#    c(
#      "Track type=GeneRegion",
#      "name=SGD Features",
#      "featureTypeAnnotation=feature_type",
#      "colorBy=feature_type",
#      "ORF=red",
#      "CDS=darkred",
#      "ARS=blue",
#      "default=gray"
#    ),
#    style_path
#  )
#  cat(sprintf("Gviz style template created: %s\n", style_path))
#} else {
#  cat("Skipping combined file creation\n")
#}
