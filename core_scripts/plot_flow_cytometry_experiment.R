#!/usr/bin/env Rscript
# Bootstrap phase
# Also loads OVERRIDE_PRESETS
FUNCTION_FILENAMES <- c("logging", "script_control", "file_operations")
for (function_filename in FUNCTION_FILENAMES) {
    function_filepath <- sprintf("~/lab_utils/core_scripts/functions_for_%s.R", function_filename)
    normalized_path <- normalizePath(function_filepath)
    if (!file.exists(normalized_path)) {
        stop(sprintf("[FATAL] File with functions not found: %s", normalized_path))
    }
    source(normalized_path)
}
message("Bootstrap phase completed...")
if(interactive()) {
  message("Interactive job... sourcing configuration file.")
  script_configuration_path <- "~/lab_utils/core_scripts/configuration_script_flow_cytometry.R"
  stopifnot(
    "Script configuration file does not exist. Please copy the template." =
    file.exists(script_configuration_path)
  )
  source(script_configuration_path)
  message("Configuration file sourced...")
}

# Ensure the variables expected in the script were //
# defined in the configuration file. //
# See template_script_configuration.R or script_configuration.R
required_configuration_variables <- c(
  "EXPERIMENT_ID",
  "SERIES_DIRECTORY",
  "DIRECTORY_ID",
  "OUTPUT_FORMAT",
  "OUTPUT_EXTENSION",
  "ACCEPT_CONFIGURATION",
  "SKIP_PACKAGE_CHECKS"
)
missing_variables <- required_configuration_variables[!sapply(required_configuration_variables, exists)]
if (length(missing_variables) > 0 ) {
  stop("Missing variable. Please define in 'script_configuration.R' file.",
       paste(missing_variables, collapse = ", "))
}
message("All variables defined in the configuration file...")

#---------------------------------------
# Confirm packages
#---------------------------------------
# Check for required packages -------------
REQUIRED_PACKAGES <- list(
  "dplyr" = "1.1.4",
  "flowCore" = "2.17.1",
  "ggcyto" = "1.26.4",
  "gtools" = "3.9.5",
  "BH" = "1.87.0.1",
  "ggplot2" = "3.3.0",
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
message("All required packages available...")
#-------------------------------------------------------------------------------
# Setup directories
#-------------------------------------------------------------------------------
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

# WARNING: I think I need to subset this if I have more than //
# one experiment in the series directory.
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
    pattern = paste0(EXPERIMENT_ID, "\\_flow_cytometry_config.R$"),
    recursive = FALSE,
    full.names = TRUE,
    include.dirs = FALSE
)

stopifnot(
  "Only one xit file expected." =
    length(XIT_FILEPATH) == 1,
  "Only one metadata_df file expected." =
    length(METADATA_FILE_PATH) == 1,
  "Only one experiment directory expected." =
    length(EXPERIMENT_DIRECTORY_PATH) == 1,
  "Only one configuration file expected." =
    length(CONFIG_FILE_PATH) == 1
)
message("Directories found...")

# Load *_CONFIG variables ---------
source(CONFIG_FILE_PATH)
required_configs <- c("EXPERIMENT_CONFIG", "RUNTIME_CONFIG")
validate_configs(required_configs)
invisible(lapply(required_configs, function(config) {
    print_config_settings(get(config), title = config)
}))

EXPECTED_SAMPLES <- EXPERIMENT_CONFIG$EXPECTED_SAMPLES
stopifnot(
  "Script EXPERIMENT_ID is not the same as CONFIG EXPERIMENT_ID" =
    EXPERIMENT_ID == EXPERIMENT_CONFIG$METADATA_df$EXPERIMENT_ID,
  "Only one xit file expected." =
    length(XIT_FILEPATH) == 1,
  "Number of fcs files should be the same as EXPECTED_SAMPLES." =
    length(FCS_FILE_PATHS) == EXPECTED_SAMPLES
)

SUBDIRS <- c("processed_data", "plots")
OUTPUT_DIRS <- sapply(SUBDIRS, function(SUBDIR) {
    subdirectory <- file.path(EXPERIMENT_DIRECTORY_PATH, SUBDIR)
    dir.create(subdirectory, recursive = TRUE, showWarnings = FALSE)
    return(subdirectory)
})
PLOT_OUTPUT_DIR <- OUTPUT_DIRS[["plots"]]

message("Variables and directories initialized...")
#-------------------------------------------------------------------------------
# Load sample metadata_df
#-------------------------------------------------------------------------------
metadata_df <- read.csv(METADATA_FILE_PATH, stringsAsFactors = FALSE)
#metadata_categories <- EXPERIMENT_CONFIG$CATEGORIES
#metadata_df$experiment_id <- EXPERIMENT_CONFIG$METADATA_df$EXPERIMENT_ID
metadata_df$file_paths <- gtools::mixedsort(FCS_FILE_PATHS)
# Convert the columns of the metadata_df to factors.
for (col_name in intersect(names(EXPERIMENT_CONFIG$CATEGORIES), colnames(metadata_df))) {
  metadata_df[[col_name]] <- factor(
    metadata_df[[col_name]],
    levels = EXPERIMENT_CONFIG$CATEGORIES[[col_name]],
    ordered = TRUE
  )
}

stopifnot(
  "Metadata_df does not have expected number of rows." =
    nrow(metadata_df) == EXPECTED_SAMPLES
)

message("Finished metadata_df processing...")
#---------------------------------------
# MAIN
#---------------------------------------
library(dplyr)
library(flowCore)
library(ggcyto)
library(gtools)
# Setup for processing --------
#DEPENDENT_COLUMN <- EXPERIMENT_CONFIG$FACET_FACTOR

# Set in the configuration with CONTROL_COLUMNS variable!
if (!is.null(EXPERIMENT_CONFIG$CONTROL_COLUMNS)) {
  CONTROL_COLUMNS <- EXPERIMENT_CONFIG$CONTROL_COLUMNS
} else {
  CONTROL_COLUMNS <- setdiff(
    names(EXPERIMENT_CONFIG$CATEGORIES),
    "timepoints"
    #EXPERIMENT_CONFIG$FACET_FACTOR
  )
}

CHANNELS_TO_PLOT <- c("FL1-A", "FSC-A", "SSC-A")

flow_set <- flowCore::read.flowSet(files = metadata_df$file_paths)

stopifnot(
  "Number of rows in metadata_df does not equal number of flow experiment." =
    length(flow_set) == nrow(metadata_df)
)
message("Loaded all flow cytometry files...")

# Assign metadata_df
# metadata has to have the same name as the flow set
rownames(metadata_df) <- basename(metadata_df$file_paths)
flowCore::pData(flow_set) <- metadata_df

# Filter and save all datasets ----------
# filtered_flow_set is a list, not an S4 class like flow_set
# Convert if needed: as(filtered_flow_frames, "flowSet")
# Apply IQR-based outlier filtering to each sample in the flow cytometry dataset
message("Filtering flow set using IQR...")
filtered_flow_set <- lapply(seq_along(flow_set), function(sample_index) {
  flow_frame <- flow_set[[sample_index]]
  event_data <- flowCore::exprs(flow_frame)
  sample_name <- flowCore::sampleNames(flow_set)[sample_index]

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
  flowCore::exprs(flow_frame) <- event_data[events_to_keep, ]

  # Report filtering results
  cat("Sample ", sample_name, " retained ", sum(events_to_keep), "/", 
      nrow(event_data), "events (", 
      round(100*sum(events_to_keep)/nrow(event_data), 1), "%)\n", sep=""
  )
  stopifnot(
    "Filtered events exceed original count!" = 
      sum(events_to_keep) <= nrow(event_data)
  )

  return(flow_frame)
})

# Convert list back to flowSet S4 class and restore metadata_df
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
channel_global_ranges <- do.call(what = rbind,
  args = lapply(CHANNELS_TO_PLOT, function(flow_channel) {
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
  }
))

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

# Extract metadata once for efficiency
fs_pdata <- pData(filtered_flow_set)

# Configuration: Define selection criteria as an expression
# Example: 
#OVERLAY_CONDITION <- quote(
#  cell_cycle_treatment == "async" & 
#  staining == "YES" &
#  rescue_allele == "4R"
#  # Add other criteria as needed, e.g.: rescue_allele == "WT"
#)
message("Considering overlay condition...")
OVERLAY_CONDITION <- NULL
if (!is.null(OVERLAY_CONDITION)) {
  # Validate expression variables exist
  expr_vars <- all.vars(OVERLAY_CONDITION)
  missing_vars <- setdiff(expr_vars, colnames(fs_pdata))
  if(length(missing_vars) > 0) {
    stop("Overlay condition contains invalid variables: ", 
         paste(missing_vars, collapse = ", "))
  }

  # 2. Evaluate condition safely
  overlay_idx <- tryCatch(
    eval(OVERLAY_CONDITION, envir = fs_pdata),
    error = function(e) {
      warning("Invalid overlay condition: ", e$message)
      return(logical(nrow(fs_pdata)))
    }
  )
  # 3. Check for single sample match
  if(sum(overlay_idx) != 1) {
    warning(paste("Overlay condition matches", sum(overlay_idx), 
            "samples. Requires exactly 1. Skipping overlay."))
    overlay_data <- NULL
  } else {
    # Extract matched sample data
    overlay_sample <- filtered_flow_set[[which(overlay_idx)]]
    overlay_data <- data.frame(
      `FL1-A` = exprs(overlay_sample)[, "FL1-A"],
      check.names = FALSE
    )
  }
}

#message("Considering async samples...")
## Check if async category exists before plotting
## TODO: Need to adjust how to set this condition potentially. //
## I forgot to take the async sample for one set. //
#if("async" %in% fs_pdata$cell_cycle_treatment) {
#
#  # Filter async samples using base R for compatibility
#  async_idx <- which(fs_pdata$cell_cycle_treatment == "async")
#  async_samples <- filtered_flow_set[async_idx]
#
#  # Create plot only if async samples exist
#  if(length(async_samples) > 0) {
#    #async_staining_plot <- ggcyto(async_samples, aes(x = `FL1-A`, color = staining)) +
#    async_staining_plot <- ggcyto(async_samples, aes(x = `FL1-A`)) +
#      geom_density(
#        aes(y = after_stat(scaled)),
#        size = 0.3,  # Use size instead of linewidth
#        alpha = 0.5
#      ) +
#      facet_grid(. ~ rescue_allele) +
#      #scale_color_manual(
#      #  values = c("YES" = "#E64B35", "NO" = "#4DBBD5"),
#      #  name = "Staining Status"
#      #) +
#      labs(
#        title = "Staining Effect in Asynchronous Samples",
#        x = "FL1-A Intensity",
#        y = "Density"
#      ) +
#      theme_minimal() +
#      theme(
#        panel.border = element_rect(color = "black", fill = NA, size = 0.8),
#        strip.background = element_blank(),
#        legend.position = "bottom"
#      )
#  }
#}


message("Plotting fl1a_plot...")
# Build base plot
#fl1a_plot <- ggcyto(filtered_flow_set[fs_pdata$cell_cycle_treatment != "async", ], 
fl1a_plot <- ggcyto(filtered_flow_set, aes(x = `FL1-A`)) +
  geom_density(
    aes(y = after_stat(scaled)),
    fill = "#4292C6",
    color = "#2166AC",
    alpha = 0.3,
    size = 0.3
  )

# Conditionally add overlay
if(exists("overlay_data")) {
  fl1a_plot <- fl1a_plot +
    geom_density(
      data = overlay_data,
      aes(x = `FL1-A`, y = after_stat(scaled)),
      color = "#636363",
      size = 0.4,
      alpha = 0.2,
      inherit.aes = FALSE
    ) +
    labs(subtitle = paste("With control overlay:", 
      identifier(overlay_sample)))
}

# Add common elements
fl1a_plot <- fl1a_plot +
  facet_grid(timepoints ~ group, switch = "y") +
  geom_vline(
    data = timepoint_medians,
    aes(xintercept = median_FL1A),
    color = "#08306B",
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
    y = "Timepoint (minutes)",
    x = "FL1-A Intensity"
  ) +
  theme_minimal()

# Add conditional annotation
if(exists("overlay_data")) {
  fl1a_plot <- fl1a_plot +
    annotate(
      "text",
      x = max(fl1a_global_range) * 0.95,
      y = 0.95,
      label = paste("Control Overlay:", identifier(overlay_sample)),
      color = "#636363",
      size = 2.5,
      hjust = 1
    )
}

#fl1a_plot <- ggcyto(filtered_flow_set, aes(x = `FL1-A`)) +
#  geom_density(
#    aes(y = after_stat(scaled)),
#    fill = "#4292C6",
#
#    color = "#2166AC",
#    #color = "#08306B", # First color.
#    alpha = 0.3,
#    size = 0.3
#  ) +
#  facet_grid(timepoints ~ group, switch = "y") +
#  geom_vline(
#    data = timepoint_medians,
#    aes(xintercept = median_FL1A),
#    color = "#E64B35",
#    linetype = "dashed",
#    size = 0.4
#  ) +
#  scale_x_continuous(
#    breaks = fl1a_global_range,
#    labels = format(fl1a_global_range, scientific = FALSE),
#    expand = c(0.02, 0)
#  ) +
#  labs(
#    title = "FL1-A Intensity Distribution",
#    subtitle = "By Timepoint and Experimental Condition",
#    y = "Timepoint (minutes)",
#    x = "FL1-A Intensity"
#  ) +
#  theme_minimal() +
#  theme(
#    # Panel customization
#    panel.background = element_blank(),
#    panel.grid.major = element_line(color = "grey80", size = 0.2),
#    panel.grid.minor = element_blank(),
#    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
#    panel.spacing = unit(0.3, "lines"),
#
#    strip.text.y.left = element_text(angle = 0, face = "bold"),
#
#    # Facet label formatting
#    strip.background = element_blank(),
#    strip.text.x = element_text(size = 6, angle = 0, hjust = 1, margin = margin(b = 5), face = "bold"),
#
#    # Axis formatting
#    axis.title = element_text(face = "bold", size = 9),
#    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, color = "gray30"),
#    axis.line.x = element_line(color = "gray60", size = 0.3),
#    axis.ticks.x = element_line(color = "gray60", size = 0.3),
#    axis.text.y = element_blank(),
#    axis.ticks.y = element_blank(),
#
#    # Title formatting
#    plot.title = element_text(face = "bold", size = 12, margin = margin(b = 5)),
#    plot.subtitle = element_text(size = 10, color = "gray30", margin = margin(b = 10)),
#    plot.margin = margin(t = 10, r = 15, b = 10, l = 15)
#
#  )

#message("Plotting fsca_vs_ssca_plot...")
# Need to add the plots from the cold spring harbor paper
#fsca_vs_ssca_plot <- ggcyto(filtered_flow_set, aes(x = `FSC-A`, y = `SSC-A`)) +
#  # ===== DATA VISUALIZATION =====
#  # Hexbin plot with viridis color scale
#  geom_hex(bins = 100, aes(fill = after_stat(density))) +
#  scale_fill_viridis_c(option = "plasma", name = "Density") +
#
#  # ===== FACETING =====
#  # Use facet_grid instead of facet_wrap to have timepoints appear only once per row
#  facet_grid(timepoints ~ group) +
#
#  # ===== AXIS SCALING =====
#  # Show only min and max on axes
#  scale_x_continuous(
#    limits = fsca_global_range,
#    breaks = fsca_breaks,
#    labels = format(fsca_breaks, scientific = TRUE, digits = 2),
#    expand = c(0.02, 0)
#  ) +
#  scale_y_continuous(
#    limits = ssca_global_range,
#    breaks = ssca_breaks,
#    labels = format(ssca_breaks, scientific = TRUE, digits = 2),
#    expand = c(0.02, 0)
#  ) +
#
#  # ===== LABELING =====
#  labs(
#    title = "Cell Size and Granularity Distribution",
#    subtitle = "FSC-A (cell size) vs SSC-A (cell granularity)",
#    x = "FSC-A (cell size)",
#    y = "SSC-A (granularity)"
#  ) +
#
#  # ===== THEME CUSTOMIZATION =====
#  theme_minimal() +
#  theme(
#    # Panel customization
#    panel.background = element_blank(),
#    panel.grid.major = element_line(color = "grey80", size = 0.2),
#    panel.grid.minor = element_line(color = "grey90", size = 0.1),
#    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
#    panel.spacing = unit(0.3, "lines"),
#
#    # Facet label formatting - timepoints on left only
#    strip.background = element_blank(),
#    strip.text.x = element_text(size = 6, face = "bold", margin = margin(b = 5)),
#    strip.text.y.left = element_text(angle = 0, face = "bold"),
#
#    # Axis formatting
#    axis.title = element_text(face = "bold", size = 9),
#    axis.text.y = element_text(size = 6, color = "gray30"),
#    axis.text.x = element_text(size = 6, angle = 45, color = "gray30", margin = margin(t = 5)),
#    axis.line = element_line(color = "gray60", size = 0.3),
#    axis.ticks = element_line(color = "gray60", size = 0.3),
#
#    # Legend formatting
#    legend.position = "right",
#    legend.key.size = unit(0.8, "lines"),
#    legend.title = element_text(size = 8, face = "bold"),
#    legend.text = element_text(size = 7),
#
#    # Title formatting
#    plot.title = element_text(face = "bold", size = 12, margin = margin(b = 5)),
#    plot.subtitle = element_text(size = 8, color = "gray30", margin = margin(b = 10)),
#    plot.margin = margin(t = 10, r = 15, b = 10, l = 15)
#  )

message("Saving plots...")
plot_object_names <- ls(pattern = "_plot$", envir = .GlobalEnv)
for (current_plot_name in plot_object_names) {
  message("  --- Plotting variables ---")
  message(sprintf("  Plotting: %s", current_plot_name))
  current_plot_object <- get(current_plot_name, envir = .GlobalEnv)
  if (!inherits(current_plot_object, "ggplot")) {
    message("Skipping ", current_plot_name, " - not a ggplot object")
    next
  }
  plot_output_path <- file.path(
    PLOT_OUTPUT_DIR,
    paste0(EXPERIMENT_ID, "_", current_plot_name, OUTPUT_EXTENSION)
  )
  message("  Saving to: ", plot_output_path)
  if (file.exists(plot_output_path)) {
    message("File already exists. Skipping")
    next
  }
  # Save the plot using ggsave. No worries about devices
  # Comment to DRY_RUN.
  ggsave(
    filename = plot_output_path,
    plot = current_plot_object,
      width = 10,       # Specify width
      height = 8,       # Specify height
      units = "in",     # Units for dimensions
      dpi = 300         # Resolution
    )
}

message("All plots saved...")
message("All done...")
