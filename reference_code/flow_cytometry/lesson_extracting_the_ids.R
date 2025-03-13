
# Script assumes working setup_flow_cytometry_experiment.R is working and was executed on an experiment.
# The data was collected from <INSERT_INSTRUMENT> with <INSERT_SOFTWARE>.
# We directly set the values to use a particular data set as an example for working with flow cytometry data.

# Helper function -------------
print_vars <- function(variable_names, environment_variables) {
    stopifnot(
        "variable_names must be character." = is.character(variable_names),
        "variable_names must be in environment ls." = all(variable_names %in% environment_variables)
    )
  for (variable in variable_names) {
        if(!typeof(get(variable)) == "closure") {
            cat(sprintf("%-25s = \n    %s\n", variable, paste(get(variable), collapse=", ")))
        }
  }
}

# Setup and validation -------------
# Load the dropbox path since that is where the data is stored.
# This assumes the variables is in the environment.
# Set manually here or in bash or use my_config repository
DROPBOX_PATH <- Sys.getenv("DROPBOX_PATH")
FLOW_CYTOMETRY_BRIDGE_PATH <- "Lab/Experiments/flow_cytometry"
if(DROPBOX_PATH == "") {
    message("Environmental variable DROPBOX_PATH not available.")
    message("Either set with my config directory or manually in the parse_flow_cytometry_arguments.")
    stop("!!!! DROPBOX_PATH required for proper directory setting.")
}

FLOW_CYTOMETRY_DIR <- file.path(DROPBOX_PATH, FLOW_CYTOMETRY_BRIDGE_PATH)
stopifnot("FLOW_CYTOMETRY_DIR does not exist." = dir.exists(FLOW_CYTOMETRY_DIR))

SERIES_NAME <- "250303_G1_arrest_degrade_and_release"
EXPERIMENT_ID <- "Exp_20250310_1"
SERIES_DIRECTORY <- file.path(FLOW_CYTOMETRY_DIR, SERIES_NAME)

XIT_FILEPATH <- list.files(
    path = SERIES_DIRECTORY,
    pattern = paste0(EXPERIMENT_ID, "\\.xit$"),
    recursive = FALSE,
    include.dirs = FALSE
)

PATHS_IN_SERIES_DIRECTORY <- dir(
    path = SERIES_DIRECTORY,
    pattern = EXPERIMENT_ID,
    full.names = TRUE,
    recursive = FALSE
)

EXPERIMENT_DIRECTORY_PATH <- PATHS_IN_SERIES_DIRECTORY[
    utils::file_test("-d", PATHS_IN_SERIES_DIRECTORY)
]

EXPECTED_SAMPLES <- 42
FCS_FILE_PATHS <- list.files(
    path = EXPERIMENT_DIRECTORY_PATH,
    pattern = "\\.fcs$",
    recursive = FALSE,
    full.names = TRUE,
    include.dirs = FALSE
)

METADATA_FILE_PATH <- list.files(
    path = SERIES_DIRECTORY,
    pattern = paste0(EXPERIMENT_ID, "\\_sample_grid.csv$"),
    recursive = FALSE,
    full.names = TRUE,
    include.dirs = FALSE
)

stopifnot(
    "Only one xit file expected." = length(XIT_FILEPATH) == 1,
    "Only one metadata file expected." = length(METADATA_FILE_PATH) == 1,
    "Only one experiment directory expected." = length(EXPERIMENT_DIRECTORY_PATH) == 1,
    "Number of fcs files should be the same as EXPECTED_SAMPLES." = length(FCS_FILE_PATHS) == EXPECTED_SAMPLES
)

OUTPUT_DIR <- "~/data/flow_cytometry_test"
SUBDIRS <- c("processed_data", "plots")
sapply(SUBDIRS, function(SUBDIR) {
    dir.create(file.path(OUTPUT_DIR, SUBDIR), recursive = TRUE, showWarnings = FALSE)
})

message("All variables initialized...")
# Print the variable to quickly inspect.
# print_vars(ls(), ls())

# Check for required packages -------------
required_packages <- c("flowCore", "svglite")
sapply(required_packages, requireNamespace, quietly = TRUE)
message("All required packages available...")

# Load metadata -------------
# Preprocess fcs file paths: They are not sorted because of the missing leading zero.
 ids <- sub(".*-([^-]+)\\.[^.]*$", "\\1", basename(FCS_FILE_PATHS))
# My code:
#numbers <- sub("[A-Z]{1,2}", "\\1", ids)
#letters <- sub("[^A-Z]{1,2}", "\\1", ids)
#pairs <- data.frame(
#    numbers = numbers,
#    letters = letters,
#    ids = ids
#)
#pairs$padded_numbers <- sapply(pairs$numbers, function(number){
#    padding_zeroes <- rep("0", max(nchar(pairs$numbers)) - nchar(number))
#    paste0(padding_zeroes, number)
#})
#
#reordered_pairs <- pairs[order(pairs$letters, pairs$padded_numbers), ]
#grepl(paste0("-", reordered_pairs$ids[1], ".fcs"), basename(FCS_FILE_PATHS))
#pairs$file_paths <- c()
#for (idx in 1:nrow(pairs)) {
#    is_file <- grepl(paste0("-", reordered_pairs$ids[idx], ".fcs"), basename(FCS_FILE_PATHS))
#    pairs$file_paths[idx] <- FCS_FILE_PATHS[is_file]
#}
# Separate letters and numbers more reliably
  letters_part <- gsub("[0-9]", "", ids)
  numbers_part <- as.numeric(gsub("[^0-9]", "", ids))
  
  # Create data frame with all components
  file_data <- data.frame(
    original_path = file_paths,
    basename = basenames,
    id = ids,
    letter = letters_part,
    number = numbers_part,
    stringsAsFactors = FALSE
  )
  
  # Add validation
  stopifnot(
    "Failed to extract valid IDs from some filenames" = all(nchar(file_data$id) > 0),
    "Failed to extract numbers from some IDs" = !anyNA(file_data$number)
  )
  
  # Sort the data frame directly
  sorted_data <- file_data[order(file_data$letter, file_data$number), ]
