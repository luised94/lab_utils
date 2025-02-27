# Load the required packages
library(readxl)
library(here)

# Function to extract text from an Excel sheet
extract_text_from_excel <- function(file_path) {
  tryCatch({
    wb <- read_excel(file_path, sheet = "Buffers", col_names = FALSE)
    text <- paste0(apply(wb, 1, function(row) {
      paste(row, collapse = " ")
    }), collapse = "\n")
    return(text)
  }, error = function(e) {
    cat(paste("Error processing", file_path, ":", e$message, "\n"))
    return(NULL)
  })
}

# Function to recursively search for Excel files in a directory
find_excel_files <- function(directory) {
  file_paths <- list.files(directory, pattern = ".xlsx$|.xls$", full.names = TRUE, recursive = TRUE)
  return(file_paths)
}

# Sun Nov  5 22:35:24 2023 ------------------------------
# Code
# Specify the directory where you want to search for Excel files
directory <- "C:/Users/liusm/Dropbox (MIT)/Lab/Experiments/ATP hydrolysis by Orc1_Orc4/Assays"
directory_to_search <- here(directory)

# Find all Excel files in the specified directory
excel_files <- find_excel_files(directory_to_search)

# Collect text from the "buffers" sheets in each Excel file
all_buffer_text <- character(0)

for (excel_file in excel_files) {
  text <- extract_text_from_excel(excel_file)
  if (!is.null(text)) {
    all_buffer_text <- c(all_buffer_text, text)
  }
}
# ... (Previous code for collecting text from Excel files)
no_NA_text <- gsub("NA", "", all_buffer_text)

# Remove duplicates from the collected text
unique_buffer_text <- unique(strsplit(no_NA_text, "\n", fixed = TRUE)[[1]])

# Rejoin the unique lines into a single text
unique_buffer_text <- paste(unique_buffer_text, collapse = "\n")

# Do something with the unique collected text (e.g., save it to a file)
writeLines(unique_buffer_text, con = "./output/unique_output.txt")

cat("Duplicates removed. Unique text saved to 'unique_output.txt'.\n")

# Split the text into recipe blocks
recipe_blocks <- strsplit(unique_buffer_text, "\n\n")

# Example of extracting recipe titles
recipe_titles <- lapply(recipe_blocks, function(block) {
  # Extract the first line as the title
  title_line <- trimws(unlist(strsplit(block, "\n", fixed = TRUE))[1])
  return(title_line)
})

