#' Find plot files based on criteria
#' @param base_dir character Base directory to search
#' @param experiment character Optional experiment ID
#' @param timestamp character Optional timestamp
#' @param pattern character Optional file pattern
#' @param verbose logical Print processing details
find_plot_files <- function(base_dir, experiment = NULL, timestamp = NULL, 
                          pattern = NULL, verbose = FALSE) {
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

#' Display plots interactively
#' @param files character vector of file paths
#' @param config list Viewer configuration
display_plots <- function(files, config) {
    if (length(files) == 0) {
        base::message("No files to display")
        return(invisible(NULL))
    }
    
    # Initialize device
    grDevices::dev.new(
        width = config$device$width,
        height = config$device$height
    )
    
    # Process each file
    for (i in base::seq_along(files)) {
        file <- files[i]
        
        if (DEBUG_CONFIG$verbose) {
            base::message(sprintf("\nDisplaying file %d of %d:", i, length(files)))
            base::message(base::basename(file))
        }
        
        # Display plot
        plot <- grDevices::readPicture(file)
        grDevices::replayPlot(plot)
        
        # Interactive viewing control
        if (DEBUG_CONFIG$interactive) {
            user_input <- base::readline(
                prompt = "Options: [Enter] next, 'p' previous, 'q' quit: "
            )
            
            if (user_input == "q") break
            if (user_input == "p" && i > 1) {
                i <- i - 2  # Will be incremented in next loop
            }
        } else {
            base::Sys.sleep(DEBUG_CONFIG$display_time)
        }
    }
    
    grDevices::dev.off()
}
