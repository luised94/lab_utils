# ==============================================================================
# Inspect Excel file: sheets, dimensions, column names, and preview
# ==============================================================================

if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
library(readxl)

# --- Configuration -----------------------------------------------------------
filepath_to_excel <- "/mnt/c/Users/Luised94/MIT Dropbox/Luis Martinez/Lab/Experiments/Loading/2022_12_18 Loading Assays Repeats for publication/analysis/260331_aggregate-analysis_load_wt-4r-supps_350mM-KGlut.xlsx"

output_filepath <- sub("\\.xlsx$", "_inspection.txt", filepath_to_excel)

# --- Inspect and write to text file ------------------------------------------
sink(output_filepath)

cat("=== FILE PATH ===\n")
cat(filepath_to_excel, "\n\n")

sheet_names <- excel_sheets(filepath_to_excel)
cat("=== SHEET NAMES ===\n")
cat(paste(seq_along(sheet_names), sheet_names, sep = ": "), sep = "\n")
cat("\n")

for (sheet_index in seq_along(sheet_names)) {
    current_sheet_name <- sheet_names[sheet_index]
    cat("=== SHEET:", current_sheet_name, "===\n")

    current_sheet_data <- read_excel(filepath_to_excel, sheet = current_sheet_name)

    cat("Dimensions:", nrow(current_sheet_data), "rows x", ncol(current_sheet_data), "columns\n\n")

    cat("Column names and types:\n")
    for (column_index in seq_along(current_sheet_data)) {
        column_name <- names(current_sheet_data)[column_index]
        column_class <- class(current_sheet_data[[column_index]])[1]
        cat(sprintf("  [%d] %-40s %s\n", column_index, column_name, column_class))
    }
    cat("\n")

    cat("First 6 rows:\n")
    print(head(current_sheet_data))
    cat("\n")

    cat("Last 6 rows:\n")
    print(tail(current_sheet_data))
    cat("\n\n")
}

sink()

cat("Inspection written to:\n", output_filepath, "\n")
