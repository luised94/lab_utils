#' Validate find_plot_files parameters
#' @param base_dir character Base directory path
#' @param experiment character Optional experiment ID
#' @param timestamp character Optional timestamp
#' @param pattern character Optional file pattern
#' @param verbose logical Print processing details
validate_find_plot_parameters <- function(base_dir, experiment, timestamp, 
                                        pattern, verbose) {
    # Required parameters
    stopifnot(
        "base_dir must be character" = is.character(base_dir),
        "base_dir must exist" = dir.exists(base_dir),
        "verbose must be logical" = is.logical(verbose)
    )
    
    # Optional parameters
    if (!is.null(experiment)) {
        stopifnot(
            "experiment must be character" = is.character(experiment),
            "experiment must match pattern" = 
                grepl(VIEWER_CONFIG$patterns$experiment, experiment)
        )
    }
    
    if (!is.null(timestamp)) {
        stopifnot(
            "timestamp must be character" = is.character(timestamp),
            "timestamp must match format" = 
                grepl(VIEWER_CONFIG$patterns$timestamp, timestamp)
        )
    }
    
    if (!is.null(pattern)) {
        stopifnot(
            "pattern must be character" = is.character(pattern)
        )
    }
}

#' Validate display_plots parameters
#' @param files character vector of file paths
#' @param config list Viewer configuration
validate_display_parameters <- function(files, config) {
    stopifnot(
        "files must be character vector" = is.character(files),
        "config must be list" = is.list(config),
        "config must contain device settings" = 
            !is.null(config$device) && 
            !is.null(config$device$width) && 
            !is.null(config$device$height),
        "device dimensions must be numeric" = 
            is.numeric(config$device$width) && 
            is.numeric(config$device$height),
        "device dimensions must be positive" = 
            config$device$width > 0 && config$device$height > 0
    )
    
    # Validate file existence
    missing_files <- files[!file.exists(files)]
    if (length(missing_files) > 0) {
        warning("Some files do not exist: ", 
                paste(basename(missing_files), collapse = ", "))
    }
}

#' Find plot files based on criteria
#' @param base_dir character Base directory to search
#' @param experiment character Optional experiment ID
#' @param timestamp character Optional timestamp
#' @param pattern character Optional file pattern
#' @param verbose logical Print processing details
find_plot_files <- function(base_dir, experiment = NULL, timestamp = NULL, 
                          pattern = NULL, verbose = FALSE) {

    # Validate inputs
    validate_find_plot_parameters(base_dir, experiment, timestamp, pattern, verbose)

    # Start with all SVG files
    files <- list.files(
        path = base_dir,
        pattern = VIEWER_CONFIG$patterns$svg,
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
    
    return(files)
}

display_plots <- function(files, config) {
    validate_display_parameters(files, config)
    if (length(files) == 0) {
        base::message("No files to display")
        return(base::invisible(NULL))
    }
    
    i <- 1
    while (i <= length(files)) {
        file <- files[i]
        
        if (DEBUG_CONFIG$verbose) {
            base::message(sprintf("\nDisplaying file %d of %d:", i, length(files)))
            base::message(base::basename(file))
        }
        
        tryCatch({
            # Read and display SVG
            img <- magick::image_read_svg(file)
            plot(img)
            
            # Interactive viewing control
            if (DEBUG_CONFIG$interactive) {
                user_input <- base::readline(
                    prompt = "Options: [Enter] next, 'p' previous, 'q' quit: "
                )
                
                if (user_input == "q") break
                if (user_input == "p" && i > 1) {
                    i <- i - 1
                } else {
                    i <- i + 1
                }
            } else {
                base::Sys.sleep(DEBUG_CONFIG$display_time)
                i <- i + 1
            }
        }, error = function(e) {
            base::message(sprintf("Error displaying file: %s", e$message))
            i <- i + 1  # Move to next file even if current fails
        })
    }
}
