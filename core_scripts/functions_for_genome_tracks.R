#' @title Calculate Global Y-limits for Bigwig Tracks
#' @description Calculates global y-axis limits from multiple bigwig files with padding
#' @param bigwig_files character vector of bigwig file paths
#' @param genome_range GRanges object specifying genomic region
#' @param padding_fraction numeric Fraction for range padding (default: 0.1)
#' @param verbose logical Print processing information
#' @return list containing {success, data, error} where data has y_limits
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges values
calculate_track_limits <- function(bigwig_files, genome_range, 
                                 padding_fraction = 0.1, verbose = FALSE) {
    result <- tryCatch({
        # Input validation
        stopifnot(
            "bigwig_files must be character vector" = is.character(bigwig_files),
            "bigwig_files cannot be empty" = length(bigwig_files) > 0,
            "genome_range must be GRanges" = inherits(genome_range, "GRanges"),
            "padding_fraction must be numeric" = is.numeric(padding_fraction),
            "padding_fraction must be between 0 and 1" = 
                padding_fraction >= 0 && padding_fraction <= 1,
            "verbose must be logical" = is.logical(verbose)
        )
        
        if (verbose) {
            message(sprintf("\nProcessing %d bigwig files for y-limits...", 
                          length(bigwig_files)))
        }
        
        # Initialize storage for values
        all_track_values <- c()
        processed_files <- 0
        skipped_files <- 0
        
        # Process each bigwig file
        for (bigwig_file in bigwig_files) {
            if (file.exists(bigwig_file)) {
                tryCatch({
                    # Import data for specific chromosome
                    track_data <- rtracklayer::import(
                        bigwig_file,
                        which = genome_range
                    )
                    
                    if (length(track_data) > 0) {
                        values <- GenomicRanges::values(track_data)$score
                        if (length(values) > 0) {
                            all_track_values <- c(all_track_values, values)
                            processed_files <- processed_files + 1
                        }
                    }
                }, error = function(e) {
                    if (verbose) {
                        message("Skipping ", basename(bigwig_file), ": ", e$message)
                    }
                    skipped_files <- skipped_files + 1
                })
            } else {
                if (verbose) {
                    message("File not found: ", basename(bigwig_file))
                }
                skipped_files <- skipped_files + 1
            }
        }
        
        # Calculate limits if we have values
        if (length(all_track_values) > 0) {
            y_min <- min(all_track_values, na.rm = TRUE)
            y_max <- max(all_track_values, na.rm = TRUE)
            y_range <- y_max - y_min
            y_limits <- c(
                y_min - (y_range * padding_fraction),
                y_max + (y_range * padding_fraction)
            )
            
            if (verbose) {
                message(sprintf("Processed %d files successfully", processed_files))
                message(sprintf("Skipped %d files", skipped_files))
                message(sprintf("Y-limits: [%.2f, %.2f]", y_limits[1], y_limits[2]))
            }
            
            list(
                success = TRUE,
                data = y_limits,
                error = NULL
            )
        } else {
            if (verbose) {
                message("No valid track data found for y-limit calculation")
            }
            
            list(
                success = FALSE,
                data = NULL,
                error = "No valid track data found"
            )
        }
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}
#' Create placeholder track for genome visualization
#' @title Create Placeholder Genomic Track
#' @description Creates an empty DataTrack with evenly spaced zero values
#'
#' @param sampling_rate numeric Number of base pairs per data point
#' @param chromosome_width numeric Total width of chromosome in base pairs
#' @param track_color character Color for track visualization
#' @param type character Track type ('l' for line, 'h' for histogram, etc.)
#' @param chromosome_name character Chromosome identifier
#' @param placeholder_format_name character Format string for track name
#' @param format_args character vector Arguments for track name formatting
#'
#' @return Gviz DataTrack object with placeholder data
create_placeholder_track <- function(
    sampling_rate = 100,
    chromosome_width = NULL,
    track_color,
    type = "l",
    chromosome_name,
    placeholder_format_name,
    format_args,
    verbose = FALSE
) {
    if (verbose) {
        message("\nCreating Placeholder Track:")
        message(sprintf("  Chromosome: %s", chromosome_name))
        message(sprintf("  Width: %d", chromosome_width))
        message(sprintf("  Points: %d", ceiling(chromosome_width / sampling_rate)))
    }

    # Input validation
    tryCatch({
        stopifnot(
            "sampling_rate must be positive numeric" = 
                is.numeric(sampling_rate) && sampling_rate > 0,
            "chromosome_width must be positive numeric" = 
                is.numeric(chromosome_width) && chromosome_width > 0,
            "track_color must be character" = 
                is.character(track_color) && length(track_color) == 1,
            "type must be valid track type" = 
                type %in% c("l", "h", "p", "g"),
            "chromosome_name must be character" = 
                is.character(chromosome_name) && length(chromosome_name) == 1,
            "placeholder_format_name must be character" = 
                is.character(placeholder_format_name),
            "format_args must be character vector" = 
                is.character(format_args)
        )

        if (verbose) message("  Input validation successful")

        # Create placeholder track
        num_points <- ceiling(chromosome_width / sampling_rate)
        
        if (verbose) message("  Creating track name...")
        placeholder_name <- tryCatch({
            do.call(sprintf, c(list(placeholder_format_name), format_args))
        }, error = function(e) {
            stop(sprintf("Failed to create placeholder name: %s", e$message))
        })

        if (verbose) message("  Creating empty ranges...")
        empty_ranges <- tryCatch({
            GenomicRanges::GRanges(
                seqnames = chromosome_name,
                ranges = IRanges::IRanges(
                    start = seq(1, chromosome_width, length.out = num_points),
                    width = 1
                ),
                score = rep(0, num_points),
                strand = "*"
            )
        }, error = function(e) {
            stop(sprintf("Failed to create genomic ranges: %s", e$message))
        })
        
        if (verbose) message("  Creating DataTrack...")
        track <- tryCatch({
            Gviz::DataTrack(
                empty_ranges,
                name = placeholder_name,
                type = type,
                col = track_color,
                chromosome = chromosome_name
            )
        }, error = function(e) {
            stop(sprintf("Failed to create DataTrack: %s", e$message))
        })

        if (verbose) {
            message("  Placeholder track created successfully:")
            message(sprintf("    Name: %s", placeholder_name))
            message(sprintf("    Points: %d", num_points))
        }

        list(
            success = TRUE,
            data = track,
            error = NULL
        )

    }, error = function(e) {
        if (verbose) {
            message("  ERROR creating placeholder track:")
            message(sprintf("    %s", e$message))
        }
        
        list(
            success = FALSE,
            data = NULL,
            error = as.character(e)
        )
    })
}

#' Create genomic data track for sample
#' @title Create Sample Track
#' @description Creates a DataTrack from bigwig file for sample visualization
#'
#' @param bigwig_file_path character Path to bigwig file
#' @param track_format_name character Format string for track name
#' @param format_args character vector Arguments for track name formatting
#' @param track_color character Color for track visualization
#' @param track_type character Track type (default: 'l' for line)
#' @param genomic_range GRanges object for data subsetting (optional)
#'
#' @return List containing:
#'   - success: logical indicating if track was created
#'   - data: DataTrack object if successful, NULL if not
#'   - error: error message if any
create_sample_track <- function(
    bigwig_file_path,
    track_format_name,
    format_args,
    track_color,
    track_type = "l",
    genomic_range = NULL,
    verbose = TRUE
) {
    # Input validation
    stopifnot(
        "bigwig_file_path must be character" = 
            is.character(bigwig_file_path) && length(bigwig_file_path) == 1,
        "track_format_name must be character" = 
            is.character(track_format_name) && length(track_format_name) == 1,
        "format_args must be character vector" = 
            is.character(format_args),
        "track_color must be character" = 
            is.character(track_color) && length(track_color) == 1,
        "track_type must be character" = 
            is.character(track_type) && length(track_type) == 1
    )
    
    if (verbose) {
        message("\nTrack Creation Attempt:")
        message(sprintf("  Bigwig path: %s", bigwig_file_path))
        message(sprintf("  File exists: %s", file.exists(bigwig_file_path)))
        message(sprintf("  Is NA: %s", is.na(bigwig_file_path)))
    }

    # Modify file check to be more explicit
    if (is.na(bigwig_file_path)) {
        return(list(
            success = FALSE,
            data = NULL,
            error = "Bigwig file path is NA"
        ))
    }
    
    if (!file.exists(bigwig_file_path)) {
        return(list(
            success = FALSE,
            data = NULL,
            error = sprintf("Bigwig file not found: %s", bigwig_file_path)
        ))
    }

    # Try to create track
    tryCatch({
        # Import data
        track_data <- if (is.null(genomic_range)) {
            rtracklayer::import(bigwig_file_path)
        } else {
            rtracklayer::import(bigwig_file_path, which = genomic_range)
        }

        # Create track name
        track_name <- do.call(sprintf, c(list(track_format_name), format_args))

        # Create track
        track <- Gviz::DataTrack(
            track_data,
            name = track_name,
            type = track_type,
            col = track_color
        )

        list(
            success = TRUE,
            data = track,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = as.character(e)
        )
    })
}

#' Create control track for genome visualization
#' @title Create Control Track
#' @description Creates a DataTrack from control sample bigwig file
#'
#' @param bigwig_file_path character Path to bigwig file
#' @param track_format_name character Format string for track name
#' @param format_args character vector Arguments for track name formatting
#' @param track_color character Color for track visualization
#' @param track_type character Track type (default: 'l' for line)
#' @param genomic_range GRanges object for data subsetting (optional)
#' @param control_type character Type of control (default: "Input")
#'
#' @return List containing success, data, and error information
create_control_track <- function(
    bigwig_file_path,
    track_format_name,
    format_args,
    track_color,
    track_type = "l",
    genomic_range = NULL
    #control_type = "Input"
) {
    stopifnot(
        "bigwig_file_path must be character" = 
            is.character(bigwig_file_path) && length(bigwig_file_path) == 1,
        "track_format_name must be character" = 
            is.character(track_format_name) && length(track_format_name) == 1,
        "format_args must be character vector" = 
            is.character(format_args),
        "track_color must be character" = 
            is.character(track_color) && length(track_color) == 1,
        "track_type must be character" = 
            is.character(track_type) && length(track_type) == 1
        #"control_type must be character" = 
        #    is.character(control_type) && length(control_type) == 1
    )

    # Check file existence
    if (is.na(bigwig_file_path) || !file.exists(bigwig_file_path)) {
        return(list(
            success = FALSE,
            data = NULL,
            error = sprintf("Control bigwig file not found: %s", bigwig_file_path)
        ))
    }

    # Try to create track
    tryCatch({
        # Import data
        track_data <- if (is.null(genomic_range)) {
            rtracklayer::import(bigwig_file_path)
        } else {
            rtracklayer::import(bigwig_file_path, which = genomic_range)
        }

        # Create track name
        track_name <- do.call(sprintf, c(list(track_format_name), format_args))

        # Create track with control-specific settings
        track <- Gviz::DataTrack(
            track_data,
            name = track_name,
            type = track_type,
            col = track_color,
            chromosome = genomic_range$seqnames[1],
            showTitle = TRUE
        )

        list(
            success = TRUE,
            data = track,
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = sprintf("Error creating control track: %s", as.character(e))
        )
    })
}

#' Create configuration for genome track visualization
#' @title Create Track Plot Configuration
#' @description Creates a configuration list for Gviz plotTracks with sensible defaults
#'
#' @param tracks List of Gviz track objects
#' @param chromosome Character, chromosome identifier
#' @param from Numeric, start position (default: 1)
#' @param to Numeric, end position
#' @param ylim Numeric vector of length 2, y-axis limits
#' @param title Character, plot title (optional)
#' @param visualization_params List of visual parameters to override defaults
#'
#' @return List of parameters for plotTracks
#' 
#' @examples
#' create_track_plot_config(
#'   tracks = tracks,
#'   chromosome = "chrI",
#'   to = 1000,
#'   ylim = c(0, 100),
#'   visualization_params = list(fontcolor = "blue")
#' )
create_track_plot_config <- function(
    tracks,
    chromosome,
    from = 1,
    to,
    ylim,
    title = NULL,
    visualization_params = list()
) {
    # Input validation
    stopifnot(
        "tracks must be a list" = is.list(tracks),
        "tracks cannot be empty" = length(tracks) > 0,
        "chromosome must be character" = is.character(chromosome),
        "from must be numeric" = is.numeric(from),
        "to must be numeric" = is.numeric(to),
        "ylim must be numeric vector of length 2" = is.numeric(ylim) && length(ylim) == 2,
        "visualization_params must be a list" = is.list(visualization_params)
    )

    # Validate all tracks are Gviz track objects
    valid_tracks <- all(sapply(tracks, inherits, "GdObject"))
    if (!valid_tracks) {
        stop("All elements in 'tracks' must be valid Gviz track objects (inherit from GdObject)")
    }

    # Hardcoded defaults
    default_params <- list(
        fontcolor = "black",           # Track name text color
        background.title = "white",    # Track name background
        col.border.title = "#E0E0E0",  # Border around track names
        cex.main = 1,                  # Main title size
        cex.axis = 0.8,               # Axis text size
        col.axis = "black",           # Axis color
        margin = 15,                  # Outer margin
        innerMargin = 5,              # Inner margin
        showTitle = TRUE              # Show track titles
    )
    
    # Merge with provided parameters
    viz_params <- modifyList(default_params, visualization_params)
    
    # Construct plot configuration
    c(
        list(
            trackList = tracks,
            chromosome = chromosome,
            from = from,
            to = to,
            ylim = ylim,
            main = title
        ),
        viz_params
    )
}

#' Execute track plotting with optional saving
#' @title Execute Track Plot
#' @description Displays and optionally saves a genome track plot
#'
#' @param plot_config List of plotting parameters from create_track_plot_config
#' @param save_path Character, path to save plot (optional)
#' @param save_params List of saving parameters (optional)
#' @return Invisible NULL
#'
#' @examples
#' execute_track_plot(
#'   plot_config = plot_config,
#'   save_path = "plot.svg",
#'   save_params = list(width = 10, height = 8)
#' )
execute_track_plot <- function(
    plot_config,
    save_path = NULL,
    save_params = list(),
    verbose = TRUE
) {
    # Validate inputs
    stopifnot(
        "plot_config must be a list" = is.list(plot_config),
        "plot_config must contain trackList" = !is.null(plot_config$trackList),
        "save_params must be a list" = is.list(save_params)
    )

    if (verbose) {
        message("\nGraphics Device Status:")
        message(sprintf("  Current device: %d", dev.cur()))
        message(sprintf("  Interactive session: %s", interactive()))
        message("\nPlot Operation:")
        message(sprintf("  Display Plot: %s", display_plot))
        message(sprintf("  Save Plot: %s", !is.null(save_path)))
        if (!is.null(save_path)) {
            message(sprintf("  Save Path: %s", save_path))
        }
        message("\nPlot Execution:")
        message(sprintf("  Number of tracks: %d", length(plot_config$trackList)))
        message("  Track types:")
        invisible(lapply(plot_config$trackList, function(track) {
            message(sprintf("    - %s: %s", 
                          class(track)[1], 
                          track@name))
        }))
    }

    
    # Modify plot function to be more explicit
    do_plot <- function() {
        if (verbose) message("\nExecuting plot...")
        
        # Ensure basic parameters are set
        plot_params <- c(plot_config, list(
            background.title = "white",
            showTitle = TRUE,
            cex.title = 1
        ))
        
        tryCatch({
            do.call(Gviz::plotTracks, plot_params)
            if (verbose) message("  Plot execution successful")
        }, error = function(e) {
            message("  Plot execution failed:")
            message(sprintf("    %s", as.character(e)))
            stop(e)
        })
    }

    # Validate save_path if provided
    if (!is.null(save_path)) {
    }

    # Default save parameters
    default_save_params <- list(
        width = 10,
        height = 8,
        device = "svg"
    )
    
    # If saving, handle device
    if (!is.null(save_path)) {
        tryCatch({
            # Check if directory exists and is writable
            save_dir <- dirname(save_path)
            if (!dir.exists(save_dir)) {
                stop("Directory does not exist: ", save_dir)
            }
            if (file.access(save_dir, mode = 2) != 0) {
                stop("Directory is not writable: ", save_dir)
            }
            # Normalize path
            save_path <- normalizePath(save_path, mustWork = FALSE)
        }, error = function(e) {
            stop("Invalid save path: ", e$message)
        })
        # Merge save parameters
        save_params <- modifyList(default_save_params, save_params)
        
        if (verbose) message("\nSaving plot...")
        # Open device based on file extension
        
        if (grepl("\\.svg$", save_path)) {
            svg(save_path, width = save_params$width, height = save_params$height)
            do_plot()
            dev.off()
        } else if (grepl("\\.pdf$", save_path)) {
            pdf(save_path, width = save_params$width, height = save_params$height)
            do_plot()
            dev.off()
        } else {
            stop("Unsupported file format. Use .svg or .pdf")
        }
    }
    # Handle display
    if (display_plot) {
        if (verbose) message("\nDisplaying plot...")
        
        if (dev.cur() == 1) {
            if (verbose) message("  Opening new graphics device")
            X11() # or quartz() on Mac, windows() on Windows
        }
        
        do_plot()
    }
    
    invisible(NULL)
}

#' Apply standard visual properties to genome tracks
#' @param tracks List of Gviz track objects
#' @param config List of visual properties from GENOME_TRACK_CONFIG
#' @return List of tracks with standardized properties
standardize_track_properties <- function(tracks, config) {
    lapply(tracks, function(track) {
        track@showTitle <- TRUE
        track@background.title <- config$track_background
        track@fontcolor.title <- config$track_fontcolor
        track@cex.title <- config$title_dev_size
        track
    })
}
