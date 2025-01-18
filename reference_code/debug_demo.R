library(stringi)
library(crayon)

get_width <- function() {
  width <- Sys.getenv("COLUMNS")
  if (width != "") {
    as.integer(width)
  } else {
    tryCatch({
      as.integer(system("tput cols", intern = TRUE))
    }, error = function(e) {
      80  # Default width if tput is not available
    })
  }
}

create_separator <- function(char = "=", width = get_width()) {
  stri_dup(char, width)
}

print_debug_info <- function(info_list, indent_char = "  ", title_char = "#") {
  if (!is.list(info_list) || length(info_list) == 0) {
    stop("info_list must be a non-empty list")
  }
  if (any(stri_isempty(names(info_list)))) {
    stop("All elements in info_list must be named")
  }
  
  # Get width and prepare output lines
  width <- get_width()
  
  # Prepare output lines
  output_lines <- c()
  
  # Process each item in the list, excluding "title"
  for (key in names(info_list)) {
    if (key == "title") next  # Skip processing for the title
    
    value <- info_list[[key]]
    indent_level <- stri_count_fixed(key, ".")
    indent <- stri_dup(indent_char, indent_level)
    clean_key <- stri_replace_all_regex(key, "^[.]+", "")
    
    if (is.null(value)) {
      line <- sprintf("%s%s:", indent, clean_key)
    } else {
      line <- sprintf("%s%s: %s", indent, clean_key, as.character(value))
    }
    
    # Append line to output_lines without modification yet
    output_lines <- c(output_lines, line)
  }
  
  # Calculate longest line length for box drawing
  max_length <- max(nchar(output_lines))
  
  # Create top and bottom separators based on longest line
  separator_line <- create_separator("=", max_length + nchar(title_char) + 2)
  
  # Apply padding and box effect to each processed line
  output_lines <- lapply(output_lines, function(line) {
    sprintf("%s-", stri_pad_right(line, max_length + nchar(title_char) + 1))  # Pad and add '-'
  })
  
  # Combine everything into final output with box effect
  final_output <- c(separator_line,
                    sprintf("%s %s-", title_char, stri_pad_right(info_list$title, max_length)),
                    separator_line,
                    output_lines,
                    separator_line)
    # Count lines made up entirely of '='
  equal_lines_count <- sum(grepl("^=+$", final_output))
  for (i in seq_along(final_output)) {
    line <- final_output[i]
    
    if (grepl("^=+$", line)) { 
      # Line composed entirely of '='
      cat(blue(line), "\n")
    } else if (grepl("^#", line)) {
      # Title line
      cat(red(line), "\n")
    } else if (i < equal_lines_count + 1) { 
      cat(blue(line), "\n")
    } else if (i < length(final_output)) { 
      # All other lines until the last one
      cat(green(line), "\n")
    } else { 
      cat(blue(line), "\n")
    }
  }
  # Use crayon for colored output
  #cat(blue(final_output[1]), "\n")   # Top separator in blue
  #cat(yellow(final_output[2]), "\n") # Title in yellow
  #cat(blue(final_output[3]), "\n")   # Middle separator in blue
  #
  #for (line in output_lines) {
  #  cat(green(line), "\n")             # Each processed line in green
  #}
  #
  #cat(blue(final_output[length(final_output)]), "\n") # Bottom separator in blue
}

# Improved example usage
debug_info <- list(
  "title" = "Debug Information",
  "Scaling Mode" = "Some Mode",
  "Plot Generation Details" = NULL,
  ".Experiment" = "Exp001",
  ".Chromosome" = "Chr1",
  ".Sample Count" = 10,
  ".Timestamp" = Sys.time(),
  ".Normalization" = "Method1",
  "Output Configuration" = NULL,
  ".Plot Directory" = "/path/to/output",
  ".Filename" = "plot.png",
  ".Full Path" = "/path/to/output/plot.png"
)

print_debug_info(debug_info)

library(futile.logger)

# Set up the logger
#flog.appender(appender.file("logfile.txt"))  # Log to a file
flog.threshold(DEBUG)  # Set minimum log level to INFO

# Example of logging at different levels
flog.debug("This is a debug message.")
flog.info("This is an info message.")
flog.warn("This is a warning message.")
flog.error("This is an error message.")

library(futile.logger)
library(stringi)
library(crayon)

# Set up the logger to log to a file
flog.appender(appender.file("debug_log.txt"))  # Change filename as needed
flog.threshold(DEBUG)  # Set minimum log level to DEBUG

get_width <- function() {
  width <- Sys.getenv("COLUMNS")
  if (width != "") {
    as.integer(width)
  } else {
    tryCatch({
      as.integer(system("tput cols", intern = TRUE))
    }, error = function(e) {
      80  # Default width if tput is not available
    })
  }
}

create_separator <- function(char = "=", width = get_width()) {
  stri_dup(char, width)
}

print_debug_info <- function(info_list, indent_char = "  ", title_char = "#") {
  # Validate input
  if (!is.list(info_list) || length(info_list) == 0) {
    flog.error("info_list must be a non-empty list")
    stop("info_list must be a non-empty list")
  }
  
  if (any(sapply(names(info_list), nchar) == 0)) {
    flog.error("All elements in info_list must be named")
    stop("All elements in info_list must be named")
  }
  
  # Log the start of the debug information print process
  flog.info("Starting to print debug information.")

  # Get width and prepare output lines
  width <- get_width()
  
  # Prepare output lines
  output_lines <- c()
  
  # Process each item in the list, excluding "title"
  for (key in names(info_list)) {
    if (key == "title") next  # Skip processing for the title
    
    value <- info_list[[key]]
    indent_level <- stri_count_fixed(key, ".")
    indent <- stri_dup(indent_char, indent_level)
    clean_key <- stri_replace_all_regex(key, "^[.]+", "")
    
    if (is.null(value)) {
      line <- sprintf("%s%s:", indent, clean_key)
    } else {
      line <- sprintf("%s%s: %s", indent, clean_key, as.character(value))
    }
    
    output_lines <- c(output_lines, line)
  }
  
  # Calculate longest line length for box drawing
  max_length <- max(nchar(output_lines))
  
  # Create top and bottom separators based on longest line
  separator_line <- create_separator("=", max_length + nchar(title_char) + 2)
  
  # Apply padding and box effect to each processed line
  output_lines_with_padding <- lapply(output_lines, function(line) {
    sprintf("%s-", stri_pad_right(line, max_length + nchar(title_char) + 1)) 
  })
  
  # Combine everything into final output with box effect
  final_output <- c(separator_line,
                    sprintf("%s %s-", title_char, stri_pad_right(info_list$title, max_length)),
                    separator_line,
                    output_lines_with_padding,
                    separator_line)

  # Convert final_output to character vector for logging
  final_output_character <- as.character(final_output)

  # Log the final output directly to a file without using logging messages
  writeLines(final_output_character, con = "debug_output.txt") 

  ###
  # Color Output (only if ANSI is supported)
  ###
  use_ansi <- crayon::has_color()
  
  # Helper for conditionally coloring lines
  color_line <- function(line, color_fun) {
    if (use_ansi) {
      color_fun(line)
    } else {
      line
    }
  }
  # Count lines made up entirely of '=' for color coding in console
  equal_lines_count <- sum(grepl("^=+$", final_output))
  
  for (i in seq_along(final_output)) {
    line <- final_output[i]
    
    if (grepl("^=+$", line)) { 
      cat(color_line(line, blue), "\n")   # Lines composed entirely of '=' in blue
    } else if (grepl("^#", line)) {
      cat(color_line(line, red), "\n")    # Title line in red
    } else if (i < equal_lines_count + 1) { 
      cat(color_line(line, blue), "\n")   # All other lines until the last one in blue
    } else if (i < length(final_output)) { 
      cat(color_line(line, green), "\n")  # Each processed line in green
    } else { 
      cat(color_line(line, blue), "\n")   # Bottom separator in blue
    }
  }

  # Log completion
  flog.info("Finished printing debug information.")
}

# Example usage
debug_info <- list(
  "title" = "Debug Information",
  "Scaling Mode" = "Some Mode",
  "Plot Generation Details" = NULL,
  ".Experiment" = "Exp001",
  ".Chromosome" = "Chr1",
  ".Sample Count" = 10,
  ".Timestamp" = Sys.time(),
  ".Normalization" = "Method1",
  "Output Configuration" = NULL,
  ".Plot Directory" = "/path/to/output",
  ".Filename" = "plot.png",
  "..Full Path" = "/path/to/output/plot.png"
)

print_debug_info(debug_info)

print_config_settings <- function(config, title = NULL, verbose = TRUE, items_per_column = 5) {
    # Input validation
    stopifnot(
        "config must be a list" = is.list(config),
        "config cannot be empty" = length(config) > 0,
        "verbose must be logical" = is.logical(verbose) && length(verbose) == 1,
        "items_per_column must be a positive integer" = is.numeric(items_per_column) && items_per_column > 0
    )

    if (!verbose) {
        return(invisible(config))
    }

    # Format title
    if (is.null(title)) {
        title <- deparse(substitute(config))
    }
    
    # Get terminal width
    width <- get_width()
    
    # Create a normalized header with '=' padding
    header_format <- "=== %s "
    header <- sprintf(header_format, title)
    
    # Calculate the remaining width for '=' padding
    remaining_width <- max(0, width - nchar(header))
    
    # Append '=' characters to fill the line
    header <- paste0(header, paste(rep("=", remaining_width), collapse = ""))

    # Ensure the header does not exceed the terminal width
    header <- substr(header, 1, width)
    separator <- create_separator("=")
    #separator <- paste(rep("=", nchar(header)), collapse = "")
    
    cat(sprintf("\n%s\n", header))
    
    # Process and format each setting
    settings <- lapply(names(config), function(setting) {
        value <- config[[setting]]
        
        # Format different types of values
        formatted_value <- if (is.logical(value)) {
            if (value) "YES" else "NO"
        } else if (is.null(value)) {
            "NULL"
        } else if (is.list(value)) {
            "LIST"  # Could expand this for nested configs if needed
        } else if (is.vector(value) & length(value) > 0) {
            escaped_newlines <- gsub("\n", "\\\\n", value)
            paste0("[", paste(escaped_newlines, collapse = ", "), "]")
        } else {
            stop(sprintf(
                "Unsupported configuration value type for setting '%s':\n  Type: %s\n  Value: %s",
                setting,
                class(value)[1],
                deparse(value)
            ))
        }
        
        # Truncate the formatted value to a maximum length (e.g., 100 characters)
        max_length <- 100
        if (nchar(formatted_value) > max_length) {
            formatted_value <- paste0(substr(formatted_value, 1, max_length - 3), "...")
        }
        # Format the line with setting name and value
        sprintf("%s: %s", gsub("_", " ", toupper(setting)), formatted_value)
    })
    
    # Determine column layout based on terminal width
    max_item_width <- max(nchar(settings)) + 2  # Add padding between columns
    num_columns <- max(1, floor(width / max_item_width))  # Ensure at least one column
    
    # Split settings into rows and columns
    num_rows <- ceiling(length(settings) / num_columns)
    
    # Pad settings into a matrix with consistent column widths
    settings_matrix <- matrix("", nrow = num_rows, ncol = num_columns)
    
    for (i in seq_along(settings)) {
        row_index <- ((i - 1) %% num_rows) + 1
        col_index <- ((i - 1) %/% num_rows) + 1
        settings_matrix[row_index, col_index] <- settings[[i]]
    }
    
    # Normalize row lengths by padding each row to the maximum row length
    row_widths <- apply(settings_matrix, 1, function(row) sum(nchar(row)))
    max_row_width <- max(row_widths)
    
    normalized_matrix <- apply(settings_matrix, 1, function(row) {
        padded_row <- sapply(seq_along(row), function(j) {
            sprintf("%-*s", max_item_width, row[j])
        })
        paste(padded_row, collapse = " ")
    })
    
    # Print normalized rows
    cat(paste(normalized_matrix, collapse = "\n"), "\n")
    
    cat(separator, "\n")
    
    # Return invisibly
    invisible(config)
}

RUNTIME_CONFIG <- list(
    debug_enabled = TRUE,
    debug_interactive = FALSE,
    debug_verbose = TRUE,
    debug_validate = TRUE,
    process_single_file = FALSE,
    process_comparison = "comp_1108forNoneAndWT",
    process_chromosome = 10,
    process_batch = 10,
    process_samples_per_batch = 4,
    process_file_index = 1,
    output_save_plots = FALSE,
    output_dry_run = TRUE,
    output_display_time = 2
)



GENOME_TRACK_CONFIG <- list(
    use_custom_visualization = FALSE,  # Control flag

    # Display dimensions
    display_width = 10,
    display_height = 8,
    
    # Track Creation
    track_points_default = 1000,
    #track_show_title = TRUE,

    # Track defaults
    track_ylim = c(0, 1000),  # Default y-limits, adjust as needed
    track_sampling_rate = 100,  # Points per base pair for empty tracks
    
    # Track colors
    color_placeholder = "#cccccc",
    color_input = "#808080",
    
    # Track naming
    format_sample_track_name = "%s: %s",
    format_control_track_name = "%s: %s - %s",
    format_placeholder_track_name = "%s: %s - %s",
    format_suffix = "(No data)",
    format_genome_axis_track_name = "Chr %s Axis",

    # Labels
    label_always_show = "antibody",
    label_never_show = c("sample_id", "full_name", "short_name", "X__cf_genotype"),
    label_separator = "-",

    # File handling
    file_pattern = "consolidated_.*_sequence\\.fastq$",
    file_sample_id = "consolidated_([0-9]{5,6})_sequence\\.fastq",
    file_sample_id_from_bigwig = "processed_([0-9]{5,6})_sequence_to_S288C_(RPKM|CPM|BPM|RPGC)\\.bw",
    file_genome_pattern = "S288C_refgenome.fna",
    file_genome_directory = file.path(Sys.getenv("HOME"), "data", "REFGENS"),
    file_feature_directory = file.path(Sys.getenv("HOME"), "data", "feature_files"),
    file_feature_pattern = "eaton_peaks",

    # File Names
    filename_format_group_templates = "%s_%s_group%02d_chr%s_%s.svg",
    filename_format_comparison_templates = "%s_%s_%s_chr%s_%s.svg",
    title_group_template = paste(
        "%s",               # Title
        "Group: %s",   # Comparison ID
        "Chromosome %s (%d samples)", # Chr info
        "%s",               # Additional info
        "Normalization: %s", # Norm method
        sep = "\n"
    ),
    title_comparison_template = paste(
        "%s",               # Title
        "Comparison: %s",   # Comparison ID
        "Chromosome %s (%d samples)", # Chr info
        "%s",               # Additional info
        "Normalization: %s", # Norm method
        sep = "\n"
    ),
    # Development mode title
    #title_dev_mode = "development",  # Enum: "development" | "publication"
    #title_dev_style = 2,    # Bold
    ## Publication mode title
    #title_pub_template = "%s: Chr%s (%s)",
    #title_pub_size = 1,
    #title_pub_style = 2,    # Bold
    ## Title constraints
    #title_max_width = 40,
    #title_max_lines = 5,
    # Interactive mode
    #interactive_prompt = "Options: [Enter] next plot, 's' skip rest, 'q' quit: ",

    track_defaults_sample = list(
        showaxis = TRUE,
        showtitle = TRUE,
        type = "h",
        size = 1.2,
        background.title = "white",
        fontcolor.title = "black",
        col.border.title = "#e0e0e0",
        cex.title = 0.6,
        fontface = 1,
        title.width = 1.2
    ),

    track_defaults_placeholder = list(
        showaxis = TRUE,
        showtitle = TRUE,
        type = "h",
        size = 0.8,
        background.title = "white",
        background.panel = "#f5f5f5",    # light gray to indicate "empty"
        fontcolor.title = "black",
        col.border.title = "#e0e0e0",
        cex.title = 0.7,
        fontface = 1,
        title.width = 0.9,
        alpha = 0.5,
        grid = FALSE
        #ylim = c(0, 1)                  # fixed range for empty tracks
    ),
    track_defaults_control = list(
        showaxis = TRUE,
        showtitle = TRUE,
        type = "h",
        size = 0.8,
        background.title = "white",
        background.panel = "white",
        fontcolor.title = "black",
        col.border.title = "#e0e0e0",
        cex.title = 0.7,
        fontface = 1,
        title.width = 0.9
        #alpha = 0.8
    ),
    track_defaults_feature = list(
        showaxis = FALSE,
        showtitle = TRUE,
        size = 0.5,
        background.title = "white",
        background.panel = "#8b7355",
        fontcolor.title = "black",
        col.border.title = "#e0e0e0",
        cex.title = 0.7,
        fontface = 1,
        title.width = 0.9,
        fill = "#8b4513",
        col = "#8b4513"
    ),
    
    # New Plot Defaults
    plot_defaults = list(
        margin = 15,
        innerMargin = 5,
        spacing = 10,
        extend.left = 0,
        extend.right = 0,
        col.axis = "black",
        cex.axis = 0.8,
        cex.main = 0.8,
        fontface.main = 2,
        background.panel = "transparent"
    )
)

iENOME_TRACK_CONFIG <- list(
    use_custom_visualization = FALSE,  # Control flag

    # Display dimensions
    display_width = 10,
    display_height = 8,
    
    # Track Creation
    track_points_default = 1000,
    #track_show_title = TRUE,

    # Track defaults
    track_ylim = c(0, 1000),  # Default y-limits, adjust as needed
    track_sampling_rate = 100,  # Points per base pair for empty tracks
    
    # Track colors
    color_placeholder = "#cccccc",
    color_input = "#808080",
    
    # Track naming
    format_sample_track_name = "%s: %s",
    format_control_track_name = "%s: %s - %s",
    format_placeholder_track_name = "%s: %s - %s",
    format_suffix = "(No data)",
    format_genome_axis_track_name = "Chr %s Axis",

    # Labels
    label_always_show = "antibody",
    label_never_show = c("sample_id", "full_name", "short_name", "X__cf_genotype"),
    label_separator = "-",

    # File handling
    file_pattern = "consolidated_.*_sequence\\.fastq$",
    file_sample_id = "consolidated_([0-9]{5,6})_sequence\\.fastq",
    file_sample_id_from_bigwig = "processed_([0-9]{5,6})_sequence_to_S288C_(RPKM|CPM|BPM|RPGC)\\.bw",
    file_genome_pattern = "S288C_refgenome.fna",
    file_genome_directory = file.path(Sys.getenv("HOME"), "data", "REFGENS"),
    file_feature_directory = file.path(Sys.getenv("HOME"), "data", "feature_files"),
    file_feature_pattern = "eaton_peaks",

    # File Names
    filename_format_group_templates = "%s_%s_group%02d_chr%s_%s.svg",
    filename_format_comparison_templates = "%s_%s_%s_chr%s_%s.svg",
    title_group_template = paste(
        "%s",               # Title
        "Group: %s",   # Comparison ID
        "Chromosome %s (%d samples)", # Chr info
        "%s",               # Additional info
        "Normalization: %s", # Norm method
        sep = "\n"
    ),
    title_comparison_template = paste(
        "%s",               # Title
        "Comparison: %s",   # Comparison ID
        "Chromosome %s (%d samples)", # Chr info
        "%s",               # Additional info
        "Normalization: %s", # Norm method
        sep = "\n"
    ),
    # Development mode title
    #title_dev_mode = "development",  # Enum: "development" | "publication"
    #title_dev_style = 2,    # Bold
    ## Publication mode title
    #title_pub_template = "%s: Chr%s (%s)",
    #title_pub_size = 1,
    #title_pub_style = 2,    # Bold
    ## Title constraints
    #title_max_width = 40,
    #title_max_lines = 5,
    # Interactive mode
    #interactive_prompt = "Options: [Enter] next plot, 's' skip rest, 'q' quit: ",

    track_defaults_sample = list(
        showaxis = TRUE,
        showtitle = TRUE,
        type = "h",
        size = 1.2,
        background.title = "white",
        fontcolor.title = "black",
        col.border.title = "#e0e0e0",
        cex.title = 0.6,
        fontface = 1,
        title.width = 1.2
    ),

    track_defaults_placeholder = list(
        showaxis = TRUE,
        showtitle = TRUE,
        type = "h",
        size = 0.8,
        background.title = "white",
        background.panel = "#f5f5f5",    # light gray to indicate "empty"
        fontcolor.title = "black",
        col.border.title = "#e0e0e0",
        cex.title = 0.7,
        fontface = 1,
        title.width = 0.9,
        alpha = 0.5,
        grid = FALSE
        #ylim = c(0, 1)                  # fixed range for empty tracks
    ),
    track_defaults_control = list(
        showaxis = TRUE,
        showtitle = TRUE,
        type = "h",
        size = 0.8,
        background.title = "white",
        background.panel = "white",
        fontcolor.title = "black",
        col.border.title = "#e0e0e0",
        cex.title = 0.7,
        fontface = 1,
        title.width = 0.9
        #alpha = 0.8
    ),
    track_defaults_feature = list(
        showaxis = FALSE,
        showtitle = TRUE,
        size = 0.5,
        background.title = "white",
        background.panel = "#8b7355",
        fontcolor.title = "black",
        col.border.title = "#e0e0e0",
        cex.title = 0.7,
        fontface = 1,
        title.width = 0.9,
        fill = "#8b4513",
        col = "#8b4513"
    ),
    
    # New Plot Defaults
    plot_defaults = list(
        margin = 15,
        innerMargin = 5,
        spacing = 10,
        extend.left = 0,
        extend.right = 0,
        col.axis = "black",
        cex.axis = 0.8,
        cex.main = 0.8,
        fontface.main = 2,
        background.panel = "transparent"
    )
)
#print_config_settings(RUNTIME_CONFIG, title = "RUNTIME CONFIGURATION")
#print_config_settings(GENOME_TRACK_CONFIG, items_per_column = 20)
