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

cat("Reading input file...\n")
raw_data <- read.delim(
  file = input_file_path,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  na.strings = "",
  quote = "",
  col.names = paste0("V", 1:16)  # Ensure we have column names
)

# Immediate debug output
debug_info <- summarize_categories(raw_data, columns = c("V2", "V3", "V12", "V9"))
print_debug_info(debug_info, "INITIAL DATA SUMMARY")

# [Rest of original processing code...]
# [Include the chromosome conversion and BED processing code from previous answer]
# [Followed by final debug output after processing]
