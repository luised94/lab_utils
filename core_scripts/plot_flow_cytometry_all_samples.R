#!/usr/bin/env Rscript
################################################################################
# Flow cytometry Plot all samples
################################################################################
# PURPOSE:
#     Determine if experiment metadata and categories are correct and output the 
#     configuration file to use in analysis scripts.
# USAGE:
#     Run from command line as script. Update template configuration by setting 
#     categories, expected samples and experiment id.
# INPUTS:
#   - cli arguments --directory-path and --experiment-id
#   - Output from the setup_flow_cytometry_experiment.R and the data
#   - Once confirmed, add --accept-configuration
# OUTPUTS: PDF with plot file, log file with reference to repository state, processed data
# CONTROLS:
#   RUNTIME_CONFIG$output_dry_run
#   RUNTIME_CONFIG$debug_verbose
#   RUNTIME_CONFIG$debug_interactive
# DEPENDENCIES:
#   - R base packages only
#   - ~/lab_utils/core_scripts/bmc_config.R
# AUTHOR: Luis
# DATE: 2025-03-11
# VERSION: 1.0.0
################################################################################
# Bootstrap phase
function_filenames <- c("logging", "script_control", "file_operations")
for (function_filename in function_filenames) {
    function_filepath <- sprintf("~/lab_utils/core_scripts/functions_for_%s.R", function_filename)
    normalized_path <- normalizePath(function_filepath)
    if (!file.exists(normalized_path)) {
        stop(sprintf("[FATAL] File with functions not found: %s", normalized_path))
    }
    source(normalized_path)
}

message("All function files loaded...")
################################################################################
# Handle script arguments
################################################################################
# Parse arguments and validate configurations
description <- "Setup flow cytometry experiments"
args <- parse_flow_cytometry_arguments(description = description)
args_info <- list(
    title = "Script Configuration",
    "script.name" = get_script_name(),
    "script.description" = description
)

print_debug_info(modifyList(args_info, args))
SERIES_DIRECTORY <- args$directory_path
EXPERIMENT_ID <- args$experiment_id
ACCEPT_CONFIGURATION <- args$accept_configuration
#EXPERIMENT_DIR <- args$experiment_dir

# Experiment ID Validation
stopifnot(
    "Only one experiment id required for this script" = length(EXPERIMENT_ID) == 1
)

message("Arguments parsed...")
################################################################################
# Setup and validation
################################################################################
# Ensure files and directories exist, parse_flow_cytometry_arguments checks as well.
DROPBOX_PATH <- Sys.getenv("DROPBOX_PATH")
FLOW_CYTOMETRY_BRIDGE_PATH <- "Lab/Experiments/flow_cytometry"
if(DROPBOX_PATH == "") {
    message("Environmental variable DROPBOX_PATH not available.")
    message("Either set with my config directory or manually in the parse_flow_cytometry_arguments.")
    stop("!!!! DROPBOX_PATH required for proper directory setting.")
}

FLOW_CYTOMETRY_DIR <- file.path(DROPBOX_PATH, FLOW_CYTOMETRY_BRIDGE_PATH)
stopifnot(
    "FLOW_CYTOMETRY_DIR does not exist." = dir.exists(FLOW_CYTOMETRY_DIR)
)

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

CONFIG_FILE_PATH <- list.files(
    path = SERIES_DIRECTORY,
    pattern = paste0(EXPERIMENT_ID, "\\flow_cytometry_config.R$"),
    recursive = FALSE,
    full.names = TRUE,
    include.dirs = FALSE
)

stopifnot(
    "Only one xit file expected." = length(XIT_FILEPATH) == 1,
    "Only one metadata file expected." = length(METADATA_FILE_PATH) == 1,
    "Only one experiment directory expected." = length(EXPERIMENT_DIRECTORY_PATH) == 1,
    "Only one experiment directory expected." = length(CONFIG_FILE_PATH) == 1
)

success <- safe_source(CONFIG_FILE_PATH, verbose = TRUE)
required_configs <- c("EXPERIMENT_CONFIG", "RUNTIME_CONFIG")
validate_configs(required_configs)
invisible(lapply(required_configs, function(config) {
    print_config_settings(get(config), title = config)
}))

EXPECTED_SAMPLES <- EXPERIMENT_CONFIG$EXPECTED_SAMPLES
stopifnot(
    "Script EXPERIMENT_ID is not the same as CONFIG EXPERIMENT_ID" = EXPERIMENT_ID == EXPERIMENT_CONFIG$METADATA$EXPERIMENT_ID,
    "Only one xit file expected." = length(xit_file) == 1,
    "Number of fcs files should be the same as EXPECTED_SAMPLES." = length(FCS_FILE_PATHS) == EXPECTED_SAMPLES
)

SUBDIRS <- c("processed_data", "plots")
OUTPUT_DIRS <- sapply(SUBDIRS, function(SUBDIR) {
    subdirectory <- file.path(EXPERIMENT_DIRECTORY_PATH, SUBDIR)
    dir.create(subdirectory, recursive = TRUE, showWarnings = FALSE)
    return(subdirectory)
})

message("Variables and directories initialized...")
################################################################################
# Directory Setup and User Confirmation
################################################################################
# Handle configuration override (independent)
if (!is.null(args$override)) {
    structured_log_info("Starting override")
    override_result <- apply_runtime_override(
        config = RUNTIME_CONFIG,
        preset_name = args$override,
        preset_list = OVERRIDE_PRESETS
    )
    RUNTIME_CONFIG <- override_result$modified

    print_debug_info(
        modifyList(
            list(
                title = "Final Configuration",
                "override.mode" = override_result$mode
                ),
            RUNTIME_CONFIG
        )
    )
    message("Override complete...")
}

handle_configuration_checkpoint(
    accept_configuration = ACCEPT_CONFIGURATION,
    experiment_id = EXPERIMENT_ID
)

message("Configuration accepted...")
########################################
# Confirm packages
########################################
# Check for required packages -------------
# TODO: Add verification for package versions
REQUIRED_PACKAGES <- list(
    "dplyr" = "1.1.4",
    "flowCore" = "2.17.1",
    "ggcyto" = "1.26.4",
    "gtools" = "3.9.5",
    "BH" = "1.87.0.1",
    "ggplot2" = "3.5.1",
    "RProtoBufLib" = "2.13.1",
    "cytolib" = "2.19.3"
)
# Verify packages are available and versions are at least minimum version in list
lapply(names(REQUIRED_PACKAGES),
    function(package_name){
        # Check availability (stop if missing - essential)
        if (!requireNamespace(package_name, quietly = TRUE)) {
             stop(paste0("Required package '", package_name, "' is not installed. ",
                         "Please install using renv or if related to flow cytometry, see install instructions."),
                  call. = FALSE) # Stop execution if a core package is missing
        }

        # Get the currently installed version
        current_version <- as.character(packageVersion(package_name))
        # Get the known compatible version
        known_version <- REQUIRED_PACKAGES[[package_name]]

        # Compare versions - issue warning if they don't match exactly
        if (compareVersion(current_version, known_version) != 0) {
            warning_message <- paste0(
                "Package '", package_name, "' version installed (", current_version, ") ",
                "differs from the known compatible version (", known_version, "). ",
                "This *may* lead to unexpected behavior. ",
                "Consider synchronizing by running 'renv::restore()'."
            )
            warning(warning_message, call. = FALSE) # Use call. = FALSE for cleaner warning output
        }
}) # END lapply

# If the loop completes without stopping, all packages are available and meet minimum version requirements.
message("All required packages are available.")

################################################################################
# Load sample metadata
################################################################################
metadata <- load_and_process_experiment_metadata(
    metadata_path = metadata_path,
    categories = EXPERIMENT_CONFIG$CATEGORIES,
    column_order = EXPERIMENT_CONFIG$COLUMN_ORDER,
    sample_ids = NULL,
    stop_on_missing_columns = TRUE
)
stopifnot(
    "Number of rows in metadata must equal expected number of samples." = nrow(metadata) == EXPECTED_SAMPLES
)

metadata$file_paths <- gtools::mixedsort(FCS_FILE_PATHS)

# Setup for processing --------
DEPENDENT_COLUMN <- EXPERIMENT_CONFIG$FACET_FACTOR

# Set in the configuration with CONTROL_FACTORS variable!
CONTROL_COLUMNS <- EXPERIMENT_CONFIG$CONTROL_COLUMNS

CHANNELS_TO_PLOT <- c("FL1-A", "FSC-A", "SSC-A")
message("Metadata file loaded...")

########################################
# MAIN
########################################
message("Loading all flow cytometry files...")
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
# Apply IQR-based outlier filtering to each sample in the flow cytometry dataset
message("Filtering flow set using IQR...")
filtered_flow_set <- lapply(seq_along(flow_set), function(sample_index) {
  flow_frame <- flow_set[[sample_index]]
  event_data <- exprs(flow_frame)
  sample_name <- sampleNames(flow_set)[sample_index]

  # Initialize logical filter for events
  events_to_keep <- rep(TRUE, nrow(event_data))

  # For each channel, identify events within 1.5*IQR of Q1 and Q3
  # Then combine these filters using a logical AND operation
  events_to_keep <- Reduce(
    `&`, 
    lapply(c("FSC-A", "SSC-A", "FL1-A"), function(channel_name) {
      channel_values <- event_data[, channel_name]
      channel_quantiles <- quantile(channel_values, probs = c(0.25, 0.75), na.rm = TRUE)
      interquantile_range <- channel_quantiles[2] - channel_quantiles[1]
      lower_bound <- channel_quantiles[1] - 1.5 * interquantile_range
      upper_bound <- channel_quantiles[2] + 1.5 * interquantile_range

      #-- return
      (channel_values >= lower_bound) & (channel_values <= upper_bound)
    })
  )

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
    apply(timepoint_medians[, CONTROL_COLUMNS], 1, paste, collapse = "_"),
    levels = unique(apply(timepoint_medians[, CONTROL_COLUMNS], 1, paste, collapse = "_")),
    ordered = TRUE
)

# Calculate global ranges for each channel
channel_global_ranges <- do.call(rbind, 
    lapply(CHANNELS_TO_PLOT, function(flow_channel) {
        channel_ranges <- fsApply(
            filtered_flow_set, function(flow_frame) range(exprs(flow_frame)[, flow_channel])
        )
        # Implicit return
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

# Get global ranges for all channels
fl1a_global_range <- channel_global_ranges %>%
  dplyr::filter(channel == "FL1-A") %>%
  select(min_value, max_value) %>%
  as.numeric()

ssca_global_range <- channel_global_ranges %>%
  dplyr::filter(channel == "SSC-A") %>%
  select(min_value, max_value) %>%
  as.numeric()

fsca_global_range <- channel_global_ranges %>%
  dplyr::filter(channel == "FSC-A") %>%
  select(min_value, max_value) %>%
  as.numeric()

# Generate pretty breaks for axes (just min and max)
fsca_breaks <- c(min(fsca_global_range), max(fsca_global_range))
ssca_breaks <- c(min(ssca_global_range), max(ssca_global_range))

# Define the file names
fl1a_file_output <- file.path(
    OUTPUT_DIRS[2],
    paste(EXPERIMENT_ID, "fl1a_density_plot.pdf", collapse = "_")
)

fsca_vs_ssca_file_output <- file.path(
    OUTPUT_DIRS[2],
    paste(EXPERIMENT_ID, "fsca_vs_ssca_plot.pdf", collapse = "_")
)

message("Plotting fl1a and fsca vs ssca plots...")
fl1a_plot <- ggcyto(filtered_flow_set, aes(x = `FL1-A`)) +
  geom_density(
    aes(y = after_stat(scaled)),
    fill = "#4292C6",
 
    color = "#2166AC",
    #color = "#08306B", # First color.
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
  scale_x_continuous(
    breaks = fl1a_global_range,
    labels = format(fl1a_global_range, scientific = FALSE),
    expand = c(0.02, 0)
  ) +
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
   
    strip.text.y.left = element_text(angle = 0, face = "bold"),

    # Facet label formatting
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, angle = 0, hjust = 1, margin = margin(b = 5), face = "bold"),
    
    # Axis formatting
    axis.title = element_text(face = "bold", size = 9),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, color = "gray30"),
    axis.line.x = element_line(color = "gray60", size = 0.3),
    axis.ticks.x = element_line(color = "gray60", size = 0.3),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    # Title formatting
    plot.title = element_text(face = "bold", size = 12, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10, color = "gray30", margin = margin(b = 10)),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 15)

  )

fsca_vs_ssca_plot <- ggcyto(filtered_flow_set, aes(x = `FSC-A`, y = `SSC-A`)) +
  # ===== DATA VISUALIZATION =====
  # Hexbin plot with viridis color scale
  geom_hex(bins = 100, aes(fill = after_stat(density))) +
  scale_fill_viridis_c(option = "plasma", name = "Density") +
  
  # ===== FACETING =====
  # Use facet_grid instead of facet_wrap to have timepoints appear only once per row
  facet_grid(timepoints ~ group) +
  
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
  
  # ===== LABELING =====
  labs(
    title = "Cell Size and Granularity Distribution",
    subtitle = "FSC-A (cell size) vs SSC-A (cell granularity)",
    x = "FSC-A (cell size)",
    y = "SSC-A (granularity)"
  ) +
  
  # ===== THEME CUSTOMIZATION =====
  theme_minimal() +
  theme(
    # Panel customization
    panel.background = element_blank(),
    panel.grid.major = element_line(color = "grey80", size = 0.2),
    panel.grid.minor = element_line(color = "grey90", size = 0.1),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.spacing = unit(0.3, "lines"),
    
    # Facet label formatting - timepoints on left only
    strip.background = element_blank(),
    strip.text.x = element_text(size = 6, face = "bold", margin = margin(b = 5)),
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    
    # Axis formatting
    axis.title = element_text(face = "bold", size = 9),
    axis.text.y = element_text(size = 6, color = "gray30"),
    axis.text.x = element_text(size = 6, angle = 45, color = "gray30", margin = margin(t = 5)),
    axis.line = element_line(color = "gray60", size = 0.3),
    axis.ticks = element_line(color = "gray60", size = 0.3),
    
    # Legend formatting
    legend.position = "right",
    legend.key.size = unit(0.8, "lines"),
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    
    # Title formatting
    plot.title = element_text(face = "bold", size = 12, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 8, color = "gray30", margin = margin(b = 10)),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 15)
  )

message("Saving plots...")
# Save the plot using ggsave. No worries about devices
ggsave(
  filename = fl1a_file_output, 
  plot = fl1a_plot,
  width = 10,       # Specify width
  height = 8,       # Specify height
  units = "in",     # Units for dimensions
  dpi = 300         # Resolution
)

ggsave(
  filename = fsca_vs_ssca_file_output,
  plot = fsca_vs_ssca_plot,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300
)
message("All done...")
