# [1] Debugging Utilities ------------------------------------------------------
summarize_categories <- function(data, columns = NULL, max_categories = 10) {
  stopifnot(
    is.data.frame(data),
    nrow(data) > 0,
    is.null(columns) || all(columns %in% names(data))
  )
  
  if(is.null(columns)) columns <- names(data)
  result <- list()
  
  for(col in columns) {
    counts <- table(data[[col]], useNA = "ifany")
    counts <- sort(counts, decreasing = TRUE)
    if(length(counts) > max_categories) {
      counts <- counts[1:max_categories]
    }
    
    result[[col]] <- sprintf(
      "%s (%d unique): %s", 
      col,
      length(unique(data[[col]])),
      paste(sprintf("%s: %d", names(counts), counts), collapse = ", ")
    )
  }
  result
}

print_debug_info <- function(debug_info, title = "DEBUG INFO") {
  stopifnot(
    is.list(debug_info),
    length(debug_info) > 0,
    all(nzchar(names(debug_info)))
  )
  
  # Create separator based on terminal width
  terminal_width <- if(!interactive()) 80 else getOption("width")
  sep <- paste(rep("=", terminal_width), collapse = "")
  
  # Create header
  cat("\n", sep, "\n", sep = "")
  cat(paste("##", toupper(title)), "\n")
  cat(sep, "\n", sep = "")
  
  # Print formatted items
  for(name in names(debug_info)) {
    items <- strwrap(debug_info[[name]], width = terminal_width - 20, exdent = 25)
    cat(sprintf("%-18s", paste0("[[ ", name, " ]]")), items[1], "\n")
    if(length(items) > 1) {
      for(item in items[-1]) cat(paste(rep(" ", 20), collapse = ""), item, "\n")
    }
  }
  
  cat(sep, "\n\n", sep = "")
}

# [2] Main Processing Script ---------------------------------------------------
input_file_path <- "~/data/feature_files/240830_SGD_features.tab"

# Initial checks
stopifnot(
  "Input file does not exist" = file.exists(input_file_path),
  "Input file is not readable" = file.access(input_file_path, 4) == 0
)


cat("Reading input file with quote handling...\n")
tryCatch({
  raw_data <- read.delim(
    file = input_file_path,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN", " "),
    quote = ""  # Disable quote interpretation
  )
}, error = function(e) {
  cat("FILE READING FAILED:\n")
  cat("First 5 lines from readLines():\n")
  print(head(readLines(input_file_path), 5))
  stop(e)
})

# Immediate debug output
debug_info <- summarize_categories(raw_data, columns = c("V2", "V3", "V12", "V9"))
# Enhanced debug summary
debug_info <- list(
  "Column Summary" = sprintf("%d columns, %d rows", ncol(raw_data), nrow(raw_data)),
  "Feature Types (V2)" = summarize_categories(raw_data["V2"]),
  "Verification Status (V3)" = summarize_categories(raw_data["V3"]),
  "Strand Codes (V12)" = summarize_categories(raw_data["V12"]),
  "Chromosomes (V9)" = summarize_categories(raw_data["V9"])
)
print_debug_info(debug_info, "DATA DIAGNOSTICS")
#print_debug_info(debug_info, "INITIAL DATA SUMMARY")

cat(sprintf("Detected %d columns\n", ncol(raw_data)))
stopifnot(
  "Unexpected column count (expected 16)" = ncol(raw_data) >= 16
)

cat(sprintf("Successfully read %d rows with %d columns\n", 
            nrow(raw_data), ncol(raw_data)))

# Column processing with validation
raw_data$chromosome_number <- as.integer(raw_data[[9]])
raw_data$start <- as.integer(raw_data[[10]])
raw_data$end <- as.integer(raw_data[[11]])
raw_data$strand_code <- raw_data[[12]]
feature_types <- raw_data[[2]]

# Validation checks
stopifnot(
  "Missing chromosome numbers" = !any(is.na(raw_data$chromosome_number)),
  "Invalid start coordinates" = !any(is.na(raw_data$start)),
  "Invalid end coordinates" = !any(is.na(raw_data$end)),
  "Invalid strand codes" = all(raw_data$strand_code %in% c("W", "C"))
)
# Chromosome conversion to roman numerals
valid_chromosomes <- c(1:17)
invalid_chrom <- unique(raw_data$chromosome_number[!raw_data$chromosome_number %in% valid_chromosomes])
stopifnot(
  "Invalid chromosome numbers detected" = length(invalid_chrom) == 0
)
cat("Converting chromosome numbers to roman numerals...\n")
raw_data$chrom <- paste0(
  "chr", 
  as.roman(raw_data$chromosome_number)
)

cat("Processing genomic coordinates...\n")
raw_data$bed_start <- pmin(raw_data$start, raw_data$end) - 1L
raw_data$bed_end <- pmax(raw_data$start, raw_data$end)
raw_data$bed_strand <- ifelse(raw_data$strand_code == "W", "+", "-")

cat("Resolving feature names...\n")

# Create ORF name lookup based on shared identifiers (V1)
orf_name_lookup <- raw_data[
  raw_data$V2 == "ORF" & !is.na(raw_data$V4), 
  c("V1", "V4")
]
colnames(orf_name_lookup) <- c("source_id", "canonical_name")

merged_data <- merge(
  raw_data,
  orf_name_lookup,
  by.x = "V1",
  by.y = "source_id",
  all.x = TRUE
)

merged_data$resolved_name <- dplyr::coalesce(
  merged_data$V4,
  merged_data$canonical_name,
  paste0(merged_data$V1, "_", merged_data$V2)  # S000350094_CDS format
)

bed_data <- data.frame(
  chrom = merged_data$chrom,
  start = merged_data$bed_start,
  end = merged_data$bed_end,
  name = merged_data$resolved_name,
  score = 0,
  strand = merged_data$bed_strand,
  stringsAsFactors = FALSE
)


stopifnot(
  "Invalid BED coordinates (start >= end)" = all(bed_data$start < bed_data$end),
  "Missing chromosome names" = !any(is.na(bed_data$chrom)),
  "Empty feature names" = !any(is.na(bed_data$name) | bed_data$name == "")
)
cat("Checking for duplicate features...\n")
duplicate_keys <- paste(bed_data$chrom, bed_data$start, bed_data$end, bed_data$name)
if(any(duplicated(duplicate_keys))) {
  dup_count <- sum(duplicated(duplicate_keys))
  cat(sprintf("WARNING: Found %d duplicate features\n", dup_count))
  print(head(bed_data[duplicated(duplicate_keys), ]))
  cat("Consider adding feature type suffixes\n")
  
  # Add feature type disambiguation
  bed_data$name <- paste0(bed_data$name, "_", merged_data$V2)
}

cat("\nWriting output files:\n")
unique_types <- unique(feature_types)
cat(sprintf("Found %d feature types: %s\n", 
            length(unique_types), paste(unique_types, collapse = ", ")))

output_directory <- normalizePath(
  file.path("~", "data", "feature_files"),
  mustWork = FALSE
)
dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)

# Verify directory creation
stopifnot(
  "Output directory creation failed" = dir.exists(output_directory),
  "Output directory not writable" = file.access(output_directory, 2) == 0
)

sanitize_filename <- function(feature_type) {
  # Convert to lowercase and replace special characters
  clean_name <- tolower(feature_type)
  # Replace all non-alphanumeric characters with underscores
  clean_name <- gsub("[^a-z0-9_]+", "_", clean_name)
  # Collapse multiple underscores
  clean_name <- gsub("_+", "_", clean_name)
  # Trim leading/trailing underscores
  clean_name <- gsub("^_|_$", "", clean_name)
  return(clean_name)
}

output_manifest <- list()
used_names <- list()

cat("\nGenerating file manifest...\n")
unique_types <- unique(feature_types)
for (feature_type in unique_types) {
  # Generate base name
  base_name <- sanitize_filename(feature_type)
  file_name <- paste0(base_name, ".bed")
  
  # Handle duplicates
  if (file_name %in% names(used_names)) {
    dup_count <- used_names[[file_name]]
    file_name <- sprintf("%s_%d.bed", base_name, dup_count + 1)
    used_names[[base_name]] <- dup_count + 1
  } else {
    used_names[[file_name]] <- 1
  }
  
  output_path <- file.path(output_directory, file_name)
  
  output_manifest[[feature_type]] <- list(
    original_type = feature_type,
    file_name = file_name,
    output_path = output_path,
    count = sum(feature_types == feature_type)
  )
}

manifest_df <- do.call(rbind, lapply(output_manifest, as.data.frame))
rownames(manifest_df) <- NULL


# Interactive verification display
cat("\n=== FILE MANIFEST (SANITIZED NAMES) ===\n")
print(manifest_df[, c("original_type", "file_name", "count")])

cat(sprintf("\nTotal features: %d | Output directory: %s\n",
            sum(manifest_df$count), output_directory))

# Duplicate resolution report
duplicate_bases <- unique(gsub("_\\d+\\.bed$", "", manifest_df$file_name))
if (length(duplicate_bases) < nrow(manifest_df)) {
  cat("\nNOTE: Some features required numeric suffixes due to name conflicts:\n")
  print(manifest_df[
    grepl("_\\d+\\.bed$", manifest_df$file_name), 
    c("original_type", "file_name")
  ])
}

response <- tolower(readline(prompt = "\nProceed with writing files? (yes/no): "))
if (!response %in% c("y", "yes")) {
  cat("Aborting file write operation\n")
  quit(save = "no", status = 1)
}
