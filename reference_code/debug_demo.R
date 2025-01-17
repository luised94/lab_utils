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

  # Count lines made up entirely of '=' for color coding in console
  equal_lines_count <- sum(grepl("^=+$", final_output))
  
  for (i in seq_along(final_output)) {
    line <- final_output[i]
    
    if (grepl("^=+$", line)) { 
      cat(blue(line), "\n")   # Lines composed entirely of '=' in blue
    } else if (grepl("^#", line)) {
      cat(red(line), "\n")    # Title line in red
    } else if (i < equal_lines_count + 1) { 
      cat(blue(line), "\n")   # All other lines until the last one in blue
    } else if (i < length(final_output)) { 
      cat(green(line), "\n")  # Each processed line in green
    } else { 
      cat(blue(line), "\n")   # Bottom separator in blue
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
  ".Full Path" = "/path/to/output/plot.png"
)

print_debug_info(debug_info)
