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
            "experiment must match pattern" = 
                grepl(validation_patterns$experiment, experiment)
        )
    }
    
    if (!is.null(timestamp)) {
        stopifnot(
            "timestamp must be character" = is.character(timestamp),
            "timestamp must match format" = 
                grepl(validation_patterns$timestamp, timestamp)
        )
    }
    
    if (!is.null(additional_patterns)) {
        # Check if single pattern or list of patterns
        if (is.character(additional_patterns)) {
            stopifnot("pattern must be non-empty" = nchar(additional_patterns) > 0)
        } else if (is.list(additional_patterns)) {
            stopifnot(
                "all patterns must be character" = all(sapply(additional_patterns, is.character)),
                "all patterns must be non-empty" = all(sapply(additional_patterns, nchar) > 0)
            )
        } else {
            stop("additional_patterns must be character or list of characters")
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
    if (!is.null(additional_patterns)) {
        if (is.character(additional_patterns)) {
            additional_patterns <- list(additional_patterns)  # Convert single pattern to list
        }
        
        # Apply all patterns
        files <- Reduce(function(files, pattern) {
            files[base::grepl(pattern, base::basename(files))]
        }, additional_patterns, init = files)
        
        if (verbose) {
            base::message(sprintf("After additional pattern filters: %d files", length(files)))
        }
    }

    return(files)
}

#' Create command configuration structure
#' @description Define available commands and their descriptions
create_viewer_commands <- function() {
    list(
        "" = "next plot",
        "p" = "previous plot",
        "q" = "quit viewer",
        "h" = "show this help",
        "i" = "show current plot info",
        "l" = "list all plots"
    )
}

#' Display help information
#' @param commands list Command configuration
#' @param verbose logical Print additional information
display_help <- function(commands, verbose = FALSE) {
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
        "device dimensions must be numeric" = 
            is.numeric(device_config$width) && 
            is.numeric(device_config$height)
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
                } else {
                    i <- i + 1
                }
            } else {
                base::Sys.sleep(display_time)
                i <- i + 1
            }
        }, error = function(e) {
            base::message(sprintf("Error displaying file: %s", e$message))
            i <- i + 1  # Move to next file even if current fails
        })
    }
}
