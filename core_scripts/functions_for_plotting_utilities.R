#' Validate plot file search parameters
#' @param base_dir character Base directory to search
#' @param experiment character Optional experiment ID
#' @param timestamp character Optional timestamp
#' @param pattern character Optional file pattern
#' @param patterns list Required patterns for validation
#' @param verbose logical Print processing details
validate_find_plot_files_parameters <- function(base_dir, experiment, timestamp, pattern, 
                                   validation_patterns, additional_filtering_patterns, verbose) {
    # Basic type checking
    stopifnot(
        "base_dir must be character" = is.character(base_dir),
        "base_dir must exist" = dir.exists(base_dir),
        "patterns must be a list" = is.list(validation_patterns),
        "required patterns missing" = all(c("experiment", "timestamp", "svg") %in% names(validation_patterns)),
        "verbose must be logical" = is.logical(verbose)
    )
    
    # Optional parameter validation
    if (!is.null(experiment)) {
        stopifnot(
            "experiment must be character" = is.character(experiment),
            "experiment must match pattern" = grepl(validation_patterns$experiment, experiment)
        )
    }
    
    if (!is.null(timestamp)) {
        stopifnot(
            "timestamp must be character" = is.character(timestamp),
            "timestamp must match format" = grepl(validation_patterns$timestamp, timestamp)
        )
    }
    
    if (!is.null(additional_filtering_patterns)) {
        # Check if single pattern or list of patterns
        if (is.character(additional_filtering_patterns)) {
            stopifnot("pattern must be non-empty" = nchar(additional_filtering_patterns) > 0)
        } else if (is.list(additional_filtering_patterns)) {
            stopifnot(
                "all patterns must be character" = all(sapply(additional_filtering_patterns, is.character)),
                "all patterns must be non-empty" = all(sapply(additional_filtering_patterns, nchar) > 0)
            )
        } else {
            stop("additional_filtering_patterns must be character or list of characters")
        }
    }
    if (!is.null(pattern)) {
        stopifnot(
            "pattern must be character" = is.character(pattern)
        )
    }
}

#' Find plot files based on criteria
#' @param base_dir character Base directory to search
#' @param experiment character Optional experiment ID
#' @param timestamp character Optional timestamp
#' @param pattern character Optional file pattern
#' @param verbose logical Print processing details
find_plot_files <- function(base_dir, patterns, experiment = NULL, 
                          timestamp = NULL, pattern = NULL, additional_filtering_patterns = NULL, verbose = FALSE) {
    # Validate inputs
    validate_find_plot_files_parameters(
        base_dir = base_dir,
        experiment = experiment,
        timestamp = timestamp,
        pattern = pattern,
        validation_patterns = patterns,
        additional_filtering_patterns = additional_filtering_patterns,
        verbose = verbose
    )
    
    # Start with all SVG files
    files <- base::list.files(
        path = base_dir,
        pattern = patterns$svg,
        recursive = TRUE,
        full.names = TRUE
    )
    
    if (verbose) {
        base::message(sprintf("Found %d SVG files", length(files)))
    }
    
    # Filter by experiment
    if (!is.null(experiment)) {
        files <- files[base::grepl(experiment, base::basename(files))]
        if (verbose) {
            base::message(sprintf("After experiment filter: %d files", length(files)))
        }
    }
    
    # Filter by timestamp
    if (!is.null(timestamp)) {
        files <- files[base::grepl(timestamp, base::basename(files))]
        if (verbose) {
            base::message(sprintf("After timestamp filter: %d files", length(files)))
        }
    }
    
    # Filter by additional pattern
    if (!is.null(pattern)) {
        files <- files[base::grepl(pattern, base::basename(files))]
        if (verbose) {
            base::message(sprintf("After pattern filter: %d files", length(files)))
        }
    }
    # Filter by additional patterns
    if (!is.null(additional_filtering_patterns)) {
        if (is.character(additional_filtering_patterns)) {
            additional_filtering_patterns <- list(additional_filtering_patterns)  # Convert single pattern to list
        }
        
        # Apply all patterns
        files <- Reduce(function(files, pattern) {
            files[base::grepl(pattern, base::basename(files))]
        }, additional_filtering_patterns, init = files)
        
        if (verbose) {
            base::message(sprintf("After additional pattern filters: %d files", length(files)))
        }
    }

    return(files)
}

#' Create command configuration structure
#' @description Define available commands and their descriptions
create_viewer_commands <- function() {
    return(list(
        "PRESS ENTER" = "next plot",
        "p" = "previous plot",
        "q" = "quit viewer",
        "h" = "show this help",
        "i" = "show current plot info",
        "l" = "list all plots",
        "g" = "go to specific plot number"
    ))
}

#' Validate commands configuration
#' @param commands list Command configuration
validate_commands <- function(commands) {
    # Basic structure validation
    stopifnot(
        "commands must be a list" = is.list(commands),
        "commands cannot be empty" = length(commands) > 0
    )
    
    # Validate contents
    invalid_names <- !sapply(names(commands), is.character)
    invalid_desc <- !sapply(commands, is.character)
    
    if (any(invalid_names)) {
        stop("All command names must be character strings")
    }
    
    if (any(invalid_desc)) {
        stop("All command descriptions must be character strings")
    }
    
    # Validate command name format (optional)
    if (any(nchar(names(commands)) > 1)) {
        warning("Some command names are longer than one character")
    }
}

#' Display help information
#' @param commands list Command configuration
#' @param verbose logical Print additional information
display_help <- function(commands, verbose = FALSE) {

    # Validation
    validate_commands(commands)

    if (verbose) {
        base::message("\nViewer Commands:")
    }
    
    # Calculate maximum command length for alignment
    max_cmd_length <- max(nchar(names(commands)))
    
    # Create format string for aligned output
    format_string <- sprintf("  %%-%ds : %%s", max_cmd_length)
    
    # Display commands
    for (cmd in names(commands)) {
        # Handle empty command (Enter key) specially
        display_cmd <- if (cmd == "") "[Enter]" else cmd
        base::message(sprintf(format_string, display_cmd, commands[[cmd]]))
    }
}

#' Validate display parameters
#' @param files character vector of file paths
#' @param device_config list Device configuration settings
validate_display_parameters <- function(files, device_config) {
    stopifnot(
        "files must be character vector" = is.character(files),
        "device_config must be list" = is.list(device_config),
        "device dimensions must be numeric" = is.numeric(device_config$width) && is.numeric(device_config$height)
    )
    
    # Validate file existence
    missing_files <- files[!file.exists(files)]
    if (length(missing_files) > 0) {
        base::warning(
            sprintf("Missing files:\n%s", 
                    paste(basename(missing_files), collapse = "\n"))
        )
    }
}

# Updated display_plots function
display_plots <- function(files, device_config, interactive = TRUE, 
                         display_time = 2, verbose = FALSE) {
    # Validate inputs
    validate_display_parameters(files, device_config)
    
    if (length(files) == 0) {
        base::message("No files to display")
        return(base::invisible(NULL))
    }

    commands <- create_viewer_commands()

    if(interactive) {
        base::message("\nStarting interactive plot viewer")
        display_help(commands, verbose = TRUE)
    }
    
    i <- 1
    while (i <= length(files)) {
        file <- files[i]
        
        if (!file.exists(file)) {
            base::message(sprintf("Skipping missing file: %s", basename(file)))
            i <- i + 1
            next
        }
        
        if (verbose) {
            base::message(sprintf("\nDisplaying file %d of %d:", i, length(files)))
            base::message(base::basename(file))
        }
        
        tryCatch({
            # Read and display SVG
            img <- magick::image_read_svg(file)
            plot(img)
            
            # Interactive viewing control
            if (interactive) {
                user_input <- base::readline(
                    prompt = sprintf("Plot %d/%d [h for help]: ", 
                                   i, length(files))
                )

                # Process user input
                if (!user_input %in% c("", names(commands))) {
                    base::message(sprintf("Invalid command: '%s'", user_input))
                    display_help(commands, verbose = TRUE)
                    next
                }
                
                # Process user input
                if (user_input == "q") {
                    break
                } else if (user_input == "h") {
                    display_help(commands, verbose = TRUE)
                    next  # Stay on current plot
                } else if (user_input == "p" && i > 1) {
                    i <- i - 1
                } else if (user_input == "i") {
                    base::message(sprintf("\nCurrent plot: %s", basename(file)))
                    next  # Stay on current plot
                } else if (user_input == "l") {
                    base::message("\nAvailable plots:")
                    for (idx in seq_along(files)) {
                        base::message(sprintf("%3d: %s", 
                                            idx, basename(files[idx])))
                    }
                    next  # Stay on current plot
                } else if (user_input == "g") {
                    base::message("\nAvailable plots:")
                    for (idx in seq_along(files)) {
                        base::message(sprintf("%3d: %s", 
                                            idx, basename(files[idx])))
                    }
                    
                    # Get plot number from user
                    plot_num <- base::readline(
                        prompt = sprintf("Enter plot number (1-%d): ", length(files))
                    )
                    
                    # Validate and convert input
                    plot_num <- tryCatch({
                        num <- as.integer(plot_num)
                        if (is.na(num) || num < 1 || num > length(files)) {
                            base::message(sprintf("Invalid plot number. Must be between 1 and %d", 
                                                length(files)))
                            i  # Keep current position
                        } else {
                            num  # Use provided number
                        }
                    }, error = function(e) {
                        base::message("Invalid input. Please enter a number.")
                        i  # Keep current position
                    })
                    
                    i <- plot_num
                    next
                } else {
                    i <- i + 1
                }
            } else {
                # If not interactive, wait display time and continue.
                base::Sys.sleep(display_time)
                i <- i + 1
            }
        }, error = function(e) {
            base::message(sprintf("Error displaying file: %s", e$message))
            i <- i + 1  # Move to next file even if current fails
        })
    }
}

#' @title Create Plot Title Using Category Information
create_plot_title <- function(metadata, comparison_name, plot_info, 
                            label_result = NULL, mode = "development") {
    # Input validation
    stopifnot(
        "metadata must be data.frame" = is.data.frame(metadata),
        "comparison_name must be character" = is.character(comparison_name),
        "plot_info must be list" = is.list(plot_info),
        "mode must be development or publication" = mode %in% c("development", "publication"),
        "label_result must contain required components" = is.list(label_result) &&       all(c("labels", "categories") %in% names(label_result$data))
    )
    
    # Use existing label information
    distinguishing_cats <- label_result$data$categories$distinguishing
    varying_cats <- label_result$data$categories$varying
    all_cats <- label_result$data$categories$all_available
    
    # Find shared characteristics using existing categorization
    shared_columns <- setdiff(all_cats, c(distinguishing_cats, "sample_id"))
    # Get shared values
    shared_values <- sapply(shared_columns, function(col) {
        values <- unique(metadata[[col]])
        if (length(values) == 1) {
            sprintf("%s: %s", col, values)
        } else {
            NULL
        }
    })
    shared_values <- unlist(shared_values[!sapply(shared_values, is.null)])
    
    # Create title based on mode
    
    if (mode == "development") {
        # Create columns with explicit line breaks and fixed width
        col1 <- sprintf(
            "%-35s\n%-35s\n%-35s\n%-35s\n%-35s",
            sprintf("Experiment: %s", plot_info$experiment_id),
            sprintf("Comparison: %s", comparison_name),
            sprintf("Chromosome: %s", plot_info$chromosome),
            sprintf("Samples: %d", nrow(metadata)),
            sprintf("Norm: %s", plot_info$normalization)
        )
        
        # Format shared properties
        shared_props <- sapply(shared_categories, function(col) {
            values <- unique(metadata[[col]])
            if (length(values) == 1) {
                sprintf("%s: %s", col, values[1])
            }
        })
        shared_props <- shared_props[!sapply(shared_props, is.null)]
        
        col2 <- sprintf(
            "Shared Properties:\n%s",
            paste(strwrap(paste(shared_props, collapse = ", "), width = 35), 
                  collapse = "\n")
        )
        
        col3 <- sprintf(
            "%-35s\n%-35s",
            format(Sys.time(), "%Y-%m-%d %H:%M"),
            sprintf("Y-range: [%.2f, %.2f]", plot_info$y_limits[1], plot_info$y_limits[2])
        )
        
        # Combine with explicit column separators
        title <- paste(
            col1,
            "|",
            col2,
            "|",
            col3,
            sep = "  "
        )
    } else {
        title <- sprintf(
            "%s: %s\nChr %s",
            plot_info$experiment_id,
            sub("^comp_", "", comparison_name),
            plot_info$chromosome
        )
    }
    
    # Format title according to configuration
    formatted_title <- format_title_text(
        title,
        max_width = GENOME_TRACK_CONFIG$title$format$max_width,
        max_lines = GENOME_TRACK_CONFIG$title$format$max_lines
    )
    
    return(formatted_title)
}

format_title_text <- function(text, max_width = 40, max_lines = NULL) {
    # Input validation
    stopifnot(
        "text must be character" = is.character(text),
        "max_width must be positive number" = is.numeric(max_width) && max_width > 0,
        "max_lines must be NULL or positive number" = 
            is.null(max_lines) || (is.numeric(max_lines) && max_lines > 0)
    )
    
    # Split into lines
    lines <- strsplit(text, "\n")[[1]]
    
    # Truncate or wrap each line
    formatted_lines <- sapply(lines, function(line) {
        if (nchar(line) > max_width) {
            paste0(substr(line, 1, max_width - 3), "...")
        } else {
            line
        }
    })
    
    # Limit number of lines if specified
    if (!is.null(max_lines) && length(formatted_lines) > max_lines) {
        formatted_lines <- c(
            formatted_lines[1:(max_lines - 1)],
            "..."
        )
    }
    
    paste(formatted_lines, collapse = "\n")
}

## Create color mapping for antibodies
#unique_antibodies <- unique(sorted_metadata$antibody)
#antibody_colors <- generate_distinct_colors(length(unique_antibodies))
#names(antibody_colors) <- unique_antibodies
#
## Update GENOME_TRACK_CONFIG with dynamic colors
#GENOME_TRACK_CONFIG$track_colors <- list(
#    antibody = antibody_colors,
#    placeholder = GENOME_TRACK_CONFIG$placeholder_color  # Maintain consistent placeholder
#)

# Create color legend text
#legend_text <- sprintf(
#    "Track Colors:\n%s\n%s",
#    paste("?", names(GENOME_TRACK_CONFIG$track_colors$antibody), 
#          sprintf("(%s)", GENOME_TRACK_CONFIG$track_colors$antibody), 
#          collapse = "\n"),
#    sprintf("? No Data (%s)", GENOME_TRACK_CONFIG$placeholder_color)
#)
