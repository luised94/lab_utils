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

########################################
# Setup and validation
########################################
# Setup paths and files -------
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

# Setup constants --------
# For creating sting identifier for rapid subsetting
COLUMN_SEPARATOR <-  "\x01"
message("All variables initialized...")
# Print the variable to quickly inspect.
# print_vars(ls(), ls())

########################################
# Load packages and metadata
########################################
# Check for required packages -------------
required_packages <- c("flowCore", "svglite", "ggcyto")
sapply(required_packages, requireNamespace, quietly = TRUE)
message("All required packages available...")

# Load metadata -------------
# Use gtools to account for unpadded digits in FCS file format.
metadata <- read.csv(METADATA_FILE_PATH)
metadata$file_paths <- gtools::mixedsort(FCS_FILE_PATHS)
#print(head(metadata))

# Setup for processing --------
DEPENDENT_COLUMN <- "timepoints"

# Set in the configuration with CONTROL_FACTORS variable!
CONTROL_COLUMNS <- c("rescue_allele", "suppressor_allele", "auxin_treatment")
UNIQUE_CONTROL_COMBINATIONS <- unique(metadata[, CONTROL_COLUMNS])

# Factor safety (prevents level->integer conversion)
metadata_chr <- lapply(metadata[CONTROL_COLUMNS], as.character)
control_combos_chr <- lapply(UNIQUE_CONTROL_COMBINATIONS[CONTROL_COLUMNS], as.character)

message("Generating keys for subsetting...")
# Create composite keys for the control columns and fast vectorized matching
metadata_keys <- do.call(paste, c(metadata_chr, sep = COLUMN_SEPARATOR))
control_keys <- do.call(paste, c(control_combos_chr, sep = COLUMN_SEPARATOR))

# Validation Phase ---------------------
EXPECTED_SAMPLES_AFTER_SUBSETTING <- length(unique(metadata[, DEPENDENT_COLUMN]))
EXPECTED_SEPARATOR_COUNT <- length(CONTROL_COLUMNS) - 1
SEPARATOR_COUNTS <- lengths(gregexpr(COLUMN_SEPARATOR, metadata_keys))
MISSING_CONTROL_IDS <- control_keys[!control_keys %in% metadata_keys]
NUMBER_OF_SAMPLES_SUBSET_WITH_CONTROL_KEYS <- lengths(lapply(control_keys, function(key){ metadata[metadata_keys %in% key, ] }))

stopifnot(
    "COLUMN_SEPARATOR must not be in FCS_FILE_PATHS strings." = all(grepl(COLUMN_SEPARATOR, x = metadata$file_paths)) == FALSE,
    "Number of rows in metadata not the same length as number of paths." = nrow(metadata) == length(FCS_FILE_PATHS),
    "Number of rows in metadata not the same as expected number of samples." = nrow(metadata) == EXPECTED_SAMPLES
)

if(!all(SEPARATOR_COUNTS == EXPECTED_SEPARATOR_COUNT)) {
  invalid_keys <- which(SEPARATOR_COUNTS != EXPECTED_SEPARATOR_COUNT)
  stop(sprintf(
    "!!!! Separator count mismatch in keys. Invalid key structure in rows = \n%s !!!!", 
    toString(head(invalid_keys, 10))
  ))
}

if (length(MISSING_CONTROL_IDS) > 0) {
  stop(
    "Control group failure: ", length(MISSING_CONTROL_IDS),
    " combinations missing in metadata. First 5: \n",
    paste(head(MISSING_CONTROL_IDS, 5), collapse = ", ")
  )
}

if (!all(NUMBER_OF_SAMPLES_SUBSET_WITH_CONTROL_KEYS == EXPECTED_SAMPLES_AFTER_SUBSETTING)) {
  invalid_keys <- which(NUMBER_OF_SAMPLES_SUBSET_WITH_CONTROL_KEYS != EXPECTED_SAMPLES_AFTER_SUBSETTING)
  stop(
    "!!!! Some of the control keys do not subset the expected number of samples: ", length(MISSING_CONTROL_IDS),
    " combinations missing in metadata. First 5: \n ",
    toString(head(invalid_keys, 5), collapse = ", ")
  )
}

message("Validation for keys complete...")

########################################
# MAIN
########################################
message("Starting main script logic...")
# Initial plots --------
# Execute for each control_key
for (key_idx in 1:length(control_keys)) {
    # Setup variables ------------
    control_key <- control_keys[key_idx]
    keys_in_current_key <- metadata_keys %in% control_key
    cat(sprintf("Processing %s key\n", key_idx))
    cat(sprintf("Using key = %s \n", control_key))

    # Subset the metadata using metadata keys -----
    # Vectorized subsetting using fast %in% operator
    # Alternative: control_samples <- merge(metadata, UNIQUE_CONTROL_COMBINATIONS, by = control_cols)
    rows_to_analyze <- metadata[keys_in_current_key, ]
    cat("Dimensions of rows_to_analyze =", dim(rows_to_analyze), "\n")
    message("Sample subsetting complete...")

    # Process rows in subset metadata -----
    message("Starting row processing...")

    # Load fcs file
    subset_flowSet <- flowCore::read.flowSet(files = rows_to_analyze$file_paths)

    if (key_idx == 1) {
        message("Initially processing first control key...")
    }

}
