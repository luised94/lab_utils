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
<<<<<<< HEAD
    stopifnot(
        "variable_names must be character." = is.character(variable_names),
        "variable_names must be in environment ls." = all(variable_names %in% environment_variables)
    )
  for (variable in variable_names) {
        if(!typeof(get(variable)) == "closure") {
            cat(sprintf("%-25s = \n    %s\n", variable, paste(get(variable), collapse=", ")))
        }
=======
  stopifnot(
    "variable_names must be character." = is.character(variable_names),
    "variable_names must be in environment ls." = all(variable_names %in% environment_variables)
  )
  for (variable in variable_names) {
    if (!typeof(get(variable)) == "closure") {
      cat(sprintf("%-25s = \n    %s\n", variable, paste(get(variable), collapse = ", ")))
    }
>>>>>>> pipeline_completion
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
<<<<<<< HEAD
if(DROPBOX_PATH == "") {
    message("Environmental variable DROPBOX_PATH not available.")
    message("Either set with my config directory or manually in the parse_flow_cytometry_arguments.")
    stop("!!!! DROPBOX_PATH required for proper directory setting.")
=======
if (DROPBOX_PATH == "") {
  message("Environmental variable DROPBOX_PATH not available.")
  message("Either set with my config directory or manually in the parse_flow_cytometry_arguments.")
  stop("!!!! DROPBOX_PATH required for proper directory setting.")
>>>>>>> pipeline_completion
}

FLOW_CYTOMETRY_DIR <- file.path(DROPBOX_PATH, FLOW_CYTOMETRY_BRIDGE_PATH)
stopifnot("FLOW_CYTOMETRY_DIR does not exist." = dir.exists(FLOW_CYTOMETRY_DIR))

SERIES_NAME <- "250303_G1_arrest_degrade_and_release"
EXPERIMENT_ID <- "Exp_20250310_1"
SERIES_DIRECTORY <- file.path(FLOW_CYTOMETRY_DIR, SERIES_NAME)

XIT_FILEPATH <- list.files(
<<<<<<< HEAD
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
=======
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
>>>>>>> pipeline_completion
]

EXPECTED_SAMPLES <- 42
FCS_FILE_PATHS <- list.files(
<<<<<<< HEAD
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
=======
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
>>>>>>> pipeline_completion
)

OUTPUT_DIR <- "~/data/flow_cytometry_test"
SUBDIRS <- c("processed_data", "plots")
sapply(SUBDIRS, function(SUBDIR) {
<<<<<<< HEAD
    dir.create(file.path(OUTPUT_DIR, SUBDIR), recursive = TRUE, showWarnings = FALSE)
=======
  dir.create(file.path(OUTPUT_DIR, SUBDIR), recursive = TRUE, showWarnings = FALSE)
>>>>>>> pipeline_completion
})

# Setup constants --------
# For creating sting identifier for rapid subsetting
<<<<<<< HEAD
COLUMN_SEPARATOR <-  "\x01"
=======
COLUMN_SEPARATOR <- "\x01"
>>>>>>> pipeline_completion
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
<<<<<<< HEAD
#print(head(metadata))
=======
# print(head(metadata))
>>>>>>> pipeline_completion

# Setup for processing --------
DEPENDENT_COLUMN <- "timepoints"

# Set in the configuration with CONTROL_FACTORS variable!
CONTROL_COLUMNS <- c("rescue_allele", "suppressor_allele", "auxin_treatment")
UNIQUE_CONTROL_COMBINATIONS <- unique(metadata[, CONTROL_COLUMNS])

CHANNELS_TO_PLOT <- c("FL1-A", "FSC-A", "SSC-A")
# Category levels for proper ordering
CATEGORIES <- list(
<<<<<<< HEAD
    rescue_allele = c("NONE", "WT", "4R"),
    suppressor_allele = c("NONE", "4PS"),
    auxin_treatment = c("NO", "YES"),
    timepoints = c("0", "20", "40", "60", "80", "100", "120")
=======
  rescue_allele = c("NONE", "WT", "4R"),
  suppressor_allele = c("NONE", "4PS"),
  auxin_treatment = c("NO", "YES"),
  timepoints = c("0", "20", "40", "60", "80", "100", "120")
>>>>>>> pipeline_completion
)

########################################
# MAIN
########################################
message("Loading all flow cytometry files...")
# Load the all the flow data.
flow_set <- flowCore::read.flowSet(files = metadata$file_paths)

stopifnot(
<<<<<<< HEAD
    "Number of rows in metadata does not equal number of flow experiment. Potential missing data." = length(flow_set) == nrow(metadata)
=======
  "Number of rows in metadata does not equal number of flow experiment. Potential missing data." =
    length(flow_set) == nrow(metadata)
>>>>>>> pipeline_completion
)

# Assign metadata
rownames(metadata) <- basename(metadata$file_paths)
flowCore::pData(flow_set) <- metadata

# Convert metadata columns to ordered factors
<<<<<<< HEAD
for(col in names(CATEGORIES)) {
  if(col %in% colnames(pData(flow_set))) {
    pData(flow_set)[[col]] <- factor(pData(flow_set)[[col]],
                                    levels = CATEGORIES[[col]],
                                    ordered = TRUE)
=======
for (col in names(CATEGORIES)) {
  if (col %in% colnames(pData(flow_set))) {
    pData(flow_set)[[col]] <- factor(pData(flow_set)[[col]],
      levels = CATEGORIES[[col]],
      ordered = TRUE
    )
>>>>>>> pipeline_completion
  }
}

# Filter and save all datasets ----------
# filtered_flow_set is a list, not an S4 class like flow_set
# Convert if needed: as(filtered_flow_frames, "flowSet")
# Apply IQR-based outlier filtering to each sample in the flow cytometry dataset
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
<<<<<<< HEAD
  cat("Sample ", sample_name, " retained ", sum(events_to_keep), "/", 
      nrow(event_data), "events (", 
      round(100*sum(events_to_keep)/nrow(event_data), 1), "%)\n", sep="")
    stopifnot("Filtered events exceed original count!" = 
          sum(events_to_keep) <= nrow(event_data))
=======
  cat("Sample ", sample_name, " retained ", sum(events_to_keep), "/",
    nrow(event_data), "events (",
    round(100 * sum(events_to_keep) / nrow(event_data), 1), "%)\n",
    sep = ""
  )
  stopifnot(
    "Filtered events exceed original count!" =
      sum(events_to_keep) <= nrow(event_data)
  )
>>>>>>> pipeline_completion

  return(flow_frame)
})

# Convert list back to flowSet S4 class and restore metadata
filtered_flow_set <- as(filtered_flow_set, "flowSet")
sampleNames(filtered_flow_set) <- sampleNames(flow_set)
pData(filtered_flow_set) <- pData(flow_set)

# Calculate median values for each sample
median_values <- fsApply(filtered_flow_set, each_col, median)
median_df <- cbind(pData(flow_set), as.data.frame(median_values))

# Group by experimental factors and calculate medians for each channel
timepoint_medians <- median_df %>%
  group_by(across(all_of(c(CONTROL_COLUMNS, "timepoints")))) %>%
  summarise(median_FL1A = median(`FL1-A`))

# Add group column to timepoint_medians for proper facet matching
timepoint_medians$group <- factor(

<<<<<<< HEAD
    apply(timepoint_medians[, CONTROL_COLUMNS], 1, paste, collapse = "_"),
    levels = unique(apply(timepoint_medians[, CONTROL_COLUMNS], 1, paste, collapse = "_")),
    ordered = TRUE
=======
  apply(timepoint_medians[, CONTROL_COLUMNS], 1, paste, collapse = "_"),
  levels = unique(apply(timepoint_medians[, CONTROL_COLUMNS], 1, paste, collapse = "_")),
  ordered = TRUE
>>>>>>> pipeline_completion
)

# Calculate global ranges for each channel
channel_global_ranges <- do.call(rbind, lapply(CHANNELS_TO_PLOT, function(flow_channel) {
  channel_ranges <- fsApply(filtered_flow_set, function(flow_frame) range(exprs(flow_frame)[, flow_channel]))
  data.frame(
    channel = flow_channel,
    min_value = min(channel_ranges),
    max_value = max(channel_ranges),
    stringsAsFactors = FALSE
  )
}))

# Add a grouping column to pData as an ordered factor
pData(filtered_flow_set)$group <- factor(
  apply(pData(filtered_flow_set)[, CONTROL_COLUMNS], 1, paste, collapse = "_"),
  levels = unique(apply(pData(filtered_flow_set)[, CONTROL_COLUMNS], 1, paste, collapse = "_")),
  ordered = TRUE
)

# Get global range for FL1-A channel
fl1a_global_range <- channel_global_ranges %>%
  dplyr::filter(channel == "FL1-A") %>%
  select(min_value, max_value) %>%
  as.numeric()

# Create the plot
file_output <- "~/flow_cytometry_test/fl1a_density_plot.pdf"

fl1a_plot <- ggcyto(filtered_flow_set, aes(x = `FL1-A`)) +
  geom_density(
    aes(y = after_stat(scaled)),
    fill = "#4292C6",
<<<<<<< HEAD
 
    color = "#2166AC",
    #color = "#08306B", # First color.
=======
    color = "#2166AC",
    # color = "#08306B", # First color.
>>>>>>> pipeline_completion
    alpha = 0.3,
    size = 0.3
  ) +
  facet_grid(timepoints ~ group, switch = "y") +
  geom_vline(
    data = timepoint_medians,
    aes(xintercept = median_FL1A),
    color = "#E64B35",
    linetype = "dashed",
    size = 0.4
  ) +
<<<<<<< HEAD
scale_x_continuous(
  breaks = fl1a_global_range,
    labels = format(fl1a_global_range, scientific = FALSE),
    expand = c(0.02, 0)
) +
    ## Replace coord_cartesian() with scale_x_continuous()
    #scale_x_continuous(
    #  limits = fl1a_global_range,
    #  breaks = fl1a_global_range,
    #  labels = format(fl1a_global_range, scientific = FALSE)
    #) +
=======
  scale_x_continuous(
    breaks = fl1a_global_range,
    labels = format(fl1a_global_range, scientific = FALSE),
    expand = c(0.02, 0)
  ) +
  ## Replace coord_cartesian() with scale_x_continuous()
  # scale_x_continuous(
  #  limits = fl1a_global_range,
  #  breaks = fl1a_global_range,
  #  labels = format(fl1a_global_range, scientific = FALSE)
  # ) +
>>>>>>> pipeline_completion
  labs(
    title = "FL1-A Intensity Distribution",
    subtitle = "By Timepoint and Experimental Condition",
    y = "Timepoint (minutes)",
    x = "FL1-A Intensity"
  ) +
  theme_minimal() +
  theme(
    # Panel customization
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "grey80", size = 0.2),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.spacing = unit(0.3, "lines"),
<<<<<<< HEAD
   
=======
>>>>>>> pipeline_completion
    strip.text.y.left = element_text(angle = 0, face = "bold"),

    # Facet label formatting
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, angle = 0, hjust = 1, margin = margin(b = 5), face = "bold"),
<<<<<<< HEAD
    
=======

>>>>>>> pipeline_completion
    # Axis formatting
    axis.title = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, color = "gray30"),
    axis.line.x = element_line(color = "gray60", size = 0.3),
    axis.ticks.x = element_line(color = "gray60", size = 0.3),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
<<<<<<< HEAD
    
=======

>>>>>>> pipeline_completion
    # Title formatting
    plot.title = element_text(face = "bold", size = 12, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10, color = "gray30", margin = margin(b = 10)),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 15)
<<<<<<< HEAD

=======
>>>>>>> pipeline_completion
  )

# Save using ggsave - no need to worry about graphics devices
ggsave(
<<<<<<< HEAD
  filename = "~/flow_cytometry_test/fl1a_density_plot.pdf", 
  plot = fl1a_plot,
  width = 10,       # Specify width
  height = 8,       # Specify height
  units = "in",     # Units for dimensions
  dpi = 300         # Resolution
=======
  filename = "~/flow_cytometry_test/fl1a_density_plot.pdf",
  plot = fl1a_plot,
  width = 10, # Specify width
  height = 8, # Specify height
  units = "in", # Units for dimensions
  dpi = 300 # Resolution
>>>>>>> pipeline_completion
)

ssca_global_range <- channel_global_ranges %>%
  dplyr::filter(channel == "SSC-A") %>%
  select(min_value, max_value) %>%
  as.numeric()

fsca_global_range <- channel_global_ranges %>%
  dplyr::filter(channel == "FSC-A") %>%
  select(min_value, max_value) %>%
  as.numeric()

fsca_vs_ssca_file_output <- "~/flow_cytometry_test/fsca_vs_ssca_plot.pdf"

# Streamlined FSC-A vs SSC-A plot matching your previous aesthetic

# Generate pretty breaks for axes (just min and max)
fsca_breaks <- c(min(fsca_global_range), max(fsca_global_range))
ssca_breaks <- c(min(ssca_global_range), max(ssca_global_range))

fsca_vs_ssca_plot <- ggcyto(filtered_flow_set, aes(x = `FSC-A`, y = `SSC-A`)) +
  # ===== DATA VISUALIZATION =====
  # Hexbin plot with viridis color scale
  geom_hex(bins = 100, aes(fill = after_stat(density))) +
  scale_fill_viridis_c(option = "plasma", name = "Density") +
<<<<<<< HEAD
  
  # ===== FACETING =====
  # Use facet_grid instead of facet_wrap to have timepoints appear only once per row
  facet_grid(timepoints ~ group) +
  
=======

  # ===== FACETING =====
  # Use facet_grid instead of facet_wrap to have timepoints appear only once per row
  facet_grid(timepoints ~ group) +

>>>>>>> pipeline_completion
  # ===== AXIS SCALING =====
  # Show only min and max on axes
  scale_x_continuous(
    limits = fsca_global_range,
    breaks = fsca_breaks,
    labels = format(fsca_breaks, scientific = TRUE, digits = 2),
    expand = c(0.02, 0)
  ) +
  scale_y_continuous(
    limits = ssca_global_range,
    breaks = ssca_breaks,
    labels = format(ssca_breaks, scientific = TRUE, digits = 2),
    expand = c(0.02, 0)
  ) +
<<<<<<< HEAD
  
=======

>>>>>>> pipeline_completion
  # ===== LABELING =====
  labs(
    title = "Cell Size and Granularity Distribution",
    subtitle = "FSC-A (cell size) vs SSC-A (cell granularity)",
    x = "FSC-A (cell size)",
    y = "SSC-A (granularity)"
  ) +
<<<<<<< HEAD
  
=======

>>>>>>> pipeline_completion
  # ===== THEME CUSTOMIZATION =====
  theme_minimal() +
  theme(
    # Panel customization
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "grey80", size = 0.2),
    panel.grid.minor = element_line(color = "grey90", size = 0.1),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.spacing = unit(0.3, "lines"),
<<<<<<< HEAD
    
=======

>>>>>>> pipeline_completion
    # Facet label formatting - timepoints on left only
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, face = "bold", margin = margin(b = 5)),
    strip.text.y.left = element_text(angle = 0, face = "bold"),
<<<<<<< HEAD
    
=======

>>>>>>> pipeline_completion
    # Axis formatting
    axis.title = element_text(face = "bold", size = 9),
    axis.text.y = element_text(size = 6, color = "gray30"),
    axis.text.x = element_text(size = 6, angle = 45, color = "gray30", margin = margin(t = 5)),
    axis.line = element_line(color = "gray60", size = 0.3),
    axis.ticks = element_line(color = "gray60", size = 0.3),
<<<<<<< HEAD
    
=======

>>>>>>> pipeline_completion
    # Legend formatting
    legend.position = "right",
    legend.key.size = unit(0.8, "lines"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
<<<<<<< HEAD
    
=======

>>>>>>> pipeline_completion
    # Title formatting
    plot.title = element_text(face = "bold", size = 12, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 8, color = "gray30", margin = margin(b = 10)),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 15)
  )

# Save the plot
ggsave(
  filename = "~/flow_cytometry_test/fsca_vs_ssca_plot.pdf",
  plot = fsca_vs_ssca_plot,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300
)
