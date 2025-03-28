# Script assumes working setup_flow_cytometry_experiment.R is working and was executed on an experiment.
# The data was collected from <INSERT_INSTRUMENT> with <INSERT_SOFTWARE>.
# We directly set the values to use a particular data set as an example for working with flow cytometry data.
# Libraries used
library(dplyr)
library(flowCore)
library(ggcyto)
library(gtools)

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
required_packages <- c("flowCore", "svglite", "ggcyto", "ggplot2")
sapply(required_packages, requireNamespace, quietly = TRUE)
message("All required packages available...")
# Get current ggplot2 version
current_version <- packageVersion("ggplot2")

# Define compatible version range
min_version <- "3.3.0"
max_version <- "3.3.6"

# Check if version is compatible
if (current_version > max_version) {
  warning(paste0("Your ggplot2 version (", current_version, ") is newer than the known compatible version (", max_version, "). This may cause issues with ggcyto."))
} else if (current_version < min_version) {
  stop(paste0("Your ggplot2 version (", current_version, ") is too old. Please upgrade to at least version ", min_version))
} else {
  message(paste0("Using ggplot2 version ", current_version, " - compatible with ggcyto"))
}

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

########################################
# MAIN
########################################
message("Loading all flow cytometry files...")
# Load the all the flow data.
flow_set <- flowCore::read.flowSet(files = metadata$file_paths)

stopifnot(
    "Number of rows in metadata does not equal number of flow experiment. Potential missing data." = length(flow_set) == nrow(metadata)
)

# Assign metadata
rownames(metadata) <- basename(metadata$file_paths)
flowCore::pData(flow_set) <- metadata

# Filter and save all datasets ----------
# filtered_flow_set is a list, not an S4 class like flow_set
# Convert if needed: as(filtered_flow_frames, "flowSet")
filtered_flow_set <- lapply(seq_along(flow_set), function(sample_index) {
  flow_frame <- flow_set[[sample_index]]
  event_data <- exprs(flow_frame)
  sample_name <- sampleNames(flow_set)[sample_index]

  # Initialize logical filter for events
  events_to_keep <- rep(TRUE, nrow(event_data))

  # For each channel, identify events within 1.5*IQR of Q1 and Q3
  # Then combine these filters using a logical AND operation
  events_to_keep <- Reduce(`&`, lapply(c("FSC-A", "SSC-A", "FL1-A"), function(channel_name) {
    channel_values <- event_data[, channel_name]
    channel_quantiles <- quantile(channel_values, probs = c(0.25, 0.75), na.rm = TRUE)
    interquantile_range <- channel_quantiles[2] - channel_quantiles[1]
    lower_bound <- channel_quantiles[1] - 1.5 * interquantile_range
    upper_bound <- channel_quantiles[2] + 1.5 * interquantile_range
    (channel_values >= lower_bound) & (channel_values <= upper_bound)
  }))

  # Apply filter to retain selected events
  exprs(flow_frame) <- event_data[events_to_keep, ]

  # Report filtering results
  cat("Sample ", sample_name, " retained ", sum(events_to_keep), "/", 
      nrow(event_data), "events (", 
      round(100*sum(events_to_keep)/nrow(event_data), 1), "%)\n", sep="")
    stopifnot("Filtered events exceed original count!" = 
          sum(events_to_keep) <= nrow(event_data))

  return(flow_frame)
})

# Convert back to flow set S4 class
filtered_flow_set <- as(filtered_flow_set, "flowSet")

median_values <- fsApply(filtered_flow_set, each_col, median)
median_df <- cbind(pData(flow_set), as.data.frame(median_values))

timepoint_medians <- median_df %>%
  group_by(across(all_of(c(CONTROL_COLUMNS, "timepoints")))) %>%
  summarise(median_FL1A = median(`FL1-A`))

CHANNELS_TO_PLOT <- c("FL1-A", "FSC-A", "SSC-A")
channel_global_ranges <- do.call(rbind, lapply(CHANNELS_TO_PLOT, function(flow_channel) {
  channel_ranges <- fsApply(filtered_flow_set, function(flow_frame) range(exprs(flow_frame)[, flow_channel]))
  data.frame(
    channel = flow_channel,
    min_value = min(channel_ranges),
    max_value = max(channel_ranges),
    stringsAsFactors = FALSE
  )
}))

fl1a_global_range <- channel_global_ranges %>%
  dplyr::filter(channel == "FL1-A") %>%
  select(min_value, max_value) %>%
  as.numeric()

sampleNames(filtered_flow_set) <- sampleNames(flow_set)
pData(filtered_flow_set) <- pData(flow_set)

pData(filtered_flow_set)$group <- apply(pData(filtered_flow_set)[, CONTROL_COLUMNS], 1, 
                                       function(x) paste(x, collapse = "_"))

pData(filtered_flow_set)$group <- factor(pData(filtered_flow_set)$group)
