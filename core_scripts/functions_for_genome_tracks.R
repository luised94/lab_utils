#' @title Calculate Global Y-limits for Bigwig Tracks
#' @description Calculates global y-axis limits from multiple bigwig files with padding
#' @param bigwig_files character vector of bigwig file paths
#' @param genome_range GRanges object specifying genomic region
#' @param padding_fraction numeric Fraction for range padding (default: 0.1)
#' @param verbose logical Print processing information
#' @return list containing {success, data, error} where data has y_limits
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges values
calculate_track_limits <- function(
    bigwig_files,
    genome_range,
    padding_fraction = 0.1,
    verbose = FALSE
) {
    if (verbose) {
        message("\nCalculating Track Limits:")
        message(sprintf("  Files to process: %d", length(bigwig_files)))
        message(sprintf("  Genome range: %s:%d-%d",
                       GenomicRanges::seqnames(genome_range),
                       GenomicRanges::start(genome_range),
                       GenomicRanges::end(genome_range)))
        message(sprintf("  Padding: %.1f%%", padding_fraction * 100))
    }

    result <- tryCatch({
        # Input validation
        stopifnot(
            "bigwig_files must be character vector" = is.character(bigwig_files),
            "bigwig_files cannot be empty" = length(bigwig_files) > 0,
            "genome_range must be GRanges" = inherits(genome_range, "GRanges"),
            "padding_fraction must be numeric" = is.numeric(padding_fraction),
            "padding_fraction must be between 0 and 1" =
                padding_fraction >= 0 && padding_fraction <= 1
        )

        # Initialize storage
        all_track_values <- c()
        processed_files <- 0
        skipped_files <- 0

        # Process files
        for (bigwig_file in bigwig_files) {
            if (verbose) {
                message(sprintf("\n  Processing: %s", basename(bigwig_file)))
            }

            if (file.exists(bigwig_file)) {
                tryCatch({
                    track_data <- rtracklayer::import(
                        bigwig_file,
                        which = genome_range
                    )

                    if (length(track_data) > 0) {
                        values <- GenomicRanges::values(track_data)$score
                        if (length(values) > 0) {
                            all_track_values <- c(all_track_values, values)
                            processed_files <- processed_files + 1

                            if (verbose) {
                                message(sprintf("    Data points: %d", length(values)))
                                message(sprintf("    Range: %.2f - %.2f",
                                              min(values),
                                              max(values)))
                            }
                        }
                    }
                }, error = function(e) {
                    if (verbose) {
                        message("    ERROR: ", e$message)
                    }
                    skipped_files <- skipped_files + 1
                })
            } else {
                if (verbose) {
                    message("    File not found")
                }
                skipped_files <- skipped_files + 1
            }
        }

        # Calculate limits
        if (length(all_track_values) > 0) {
            y_min <- min(all_track_values, na.rm = TRUE)
            y_max <- max(all_track_values, na.rm = TRUE)
            y_range <- y_max - y_min
            y_limits <- c(
                max(0, y_min - (y_range * padding_fraction)),
                y_max + (y_range * padding_fraction)
            )

            if (verbose) {
                message("\nLimit Calculation Summary:")
                message(sprintf("  Files processed: %d", processed_files))
                message(sprintf("  Files skipped: %d", skipped_files))
                message(sprintf("  Total data points: %d",
                              length(all_track_values)))
                message(sprintf("  Raw range: [%.2f, %.2f]", y_min, y_max))
                message(sprintf("  Padded limits: [%.2f, %.2f]",
                              y_limits[1], y_limits[2]))
            }

            list(
                success = TRUE,
                data = y_limits,
                error = NULL
            )
        } else {
            if (verbose) {
                message("\nNo valid track data found")
                message(sprintf("  Files processed: %d", processed_files))
                message(sprintf("  Files skipped: %d", skipped_files))
            }

            list(
                success = FALSE,
                data = NULL,
                error = "No valid track data found"
            )
        }
    }, error = function(e) {
        if (verbose) {
            message("\nERROR in limit calculation:")
            message(sprintf("  %s", e$message))
        }

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
    track_type = "h",
    chromosome_name,
    placeholder_format_name,
    format_args,
    track_params = list(),
    verbose = FALSE
) {
    if (verbose) {
        message("\nCreating Placeholder Track:")
        message(sprintf("  Chromosome: %s", chromosome_name))
        message(sprintf("  Width: %d", chromosome_width))
        message(sprintf("  Points: %d", ceiling(chromosome_width / sampling_rate)))
        if (length(track_params) > 0) {
            message("  Track Parameters:")
            invisible(lapply(names(track_params), function(param) {
                message(sprintf("    %s: %s",
                              param,
                              paste(track_params[[param]], collapse = ", ")))
            }))
        }
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
            "track_type must be valid track type" =
                track_type %in% c("l", "h", "p", "g"),
            "chromosome_name must be character" =
                is.character(chromosome_name) && length(chromosome_name) == 1,
            "placeholder_format_name must be character" =
                is.character(placeholder_format_name),
            "format_args must be character vector" =
                is.character(format_args),
            "track_params must be a list" =
                is.list(track_params)
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

        # Merge default parameters with provided ones
        track_args <- c(
            list(
                empty_ranges,
                name = placeholder_name,
                type = track_type,
                col = track_color,
                chromosome = chromosome_name  # Added explicit chromosome
            ),
            track_params
        )

        if (verbose) {
            message("  Creating DataTrack with parameters:")
            message(sprintf("    Type: %s", track_type))
            message(sprintf("    Color: %s", track_color))
            message(sprintf("    Additional params: %d", length(track_params)))
        }

        track <- tryCatch({
            do.call(Gviz::DataTrack, track_args)
        }, error = function(e) {
            stop(sprintf("Failed to create DataTrack: %s", e$message))
        })


        if (verbose) {
            message("  Placeholder track created successfully:")
            message(sprintf("    Name: %s", placeholder_name))
            message(sprintf("    Points: %d", num_points))
            #message(sprintf("    Track size: %s",
                          #if(is.null(track@size)) "default" else track@size))
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
# TODO: Really have to reconsider if this needs to be a function in this manner //
# has a lot of code because of debugging, validation and return. //
create_sample_track <- function(
    bigwig_file_path,
    track_format_name,
    format_args,
    track_color,
    track_type = "l",
    genomic_range = NULL,
    track_params = list(),
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
            is.character(track_type) && length(track_type) == 1,
        "track_params must be a list" =
            is.list(track_params)
    )

    if (verbose) {
        message("\nTrack Creation Attempt:")
        message(sprintf("  Bigwig path: %s", bigwig_file_path))
        message(sprintf("  Is NA: %s", is.na(bigwig_file_path)))
        message(sprintf("  File exists: %s", file.exists(bigwig_file_path)))
        message(sprintf("  Track type: %s", track_type))
        if (length(track_params) > 0) {
            message("  Track parameters:")
            invisible(lapply(names(track_params), function(param) {
                message(sprintf("    %s: %s",
                              param,
                              paste(track_params[[param]], collapse = ", ")))
            }))
        }
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

        if (verbose) message(sprintf("  Creating track: %s", track_name))

        # Merge default parameters with provided ones
        track_args <- c(
            list(
                track_data,
                name = track_name,
                type = track_type,
                col = track_color
            ),
            track_params
        )

        # Create track with all parameters
        track <- do.call(Gviz::DataTrack, track_args)

        if (verbose) message("  Track creation successful")

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
    track_type = "h",
    genomic_range = NULL,
    track_params = list(),
    verbose = FALSE
) {
    if (verbose) {
        message("\nCreating Control Track:")
        message(sprintf("  Bigwig path: %s", bigwig_file_path))
        message(sprintf("  File exists: %s", file.exists(bigwig_file_path)))
        if (!is.null(genomic_range)) {
            message(sprintf("  Range: %s:%d-%d",
                          GenomicRanges::seqnames(genomic_range),
                          GenomicRanges::start(genomic_range),
                          GenomicRanges::end(genomic_range)))
        }
        if (length(track_params) > 0) {
            message("  Track Parameters:")
            invisible(lapply(names(track_params), function(param) {
                message(sprintf("    %s: %s",
                              param,
                              paste(track_params[[param]], collapse = ", ")))
            }))
        }
    }

    # Input validation
    tryCatch({
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
                is.character(track_type) && length(track_type) == 1,
            "track_params must be a list" =
                is.list(track_params)
        )

        if (verbose) message("  Input validation successful")

        # Check file existence
        if (is.na(bigwig_file_path) || !file.exists(bigwig_file_path)) {
            return(list(
                success = FALSE,
                data = NULL,
                error = sprintf("Control bigwig file not found: %s", bigwig_file_path)
            ))
        }

        # Import data
        if (verbose) message("  Importing bigwig data...")
        track_data <- if (is.null(genomic_range)) {
            rtracklayer::import(bigwig_file_path)
        } else {
            rtracklayer::import(bigwig_file_path, which = genomic_range)
        }

        # Create track name
        if (verbose) message("  Creating track name...")
        track_name <- tryCatch({
            do.call(sprintf, c(list(track_format_name), format_args))
        }, error = function(e) {
            stop(sprintf("Failed to create track name: %s", e$message))
        })

        # Merge default parameters with provided ones
        track_args <- c(
            list(
                track_data,
                name = track_name,
                type = track_type,
                col = track_color,
                chromosome = if (!is.null(genomic_range))
                    as.character(GenomicRanges::seqnames(genomic_range)[1])
                else NULL
            ),
            track_params
        )

        if (verbose) {
            message("  Creating DataTrack with parameters:")
            message(sprintf("    Type: %s", track_type))
            message(sprintf("    Color: %s", track_color))
            message(sprintf("    Additional params: %d", length(track_params)))
        }

        # Create track
        track <- tryCatch({
            do.call(Gviz::DataTrack, track_args)
        }, error = function(e) {
            stop(sprintf("Failed to create DataTrack: %s", e$message))
        })


        if (verbose) {
            message("  Control track created successfully:")
            message(sprintf("    Name: %s", track_name))
            message(sprintf("    Data points: %d", length(track_data)))
            #message(sprintf("    Track size: %s",
                          #if(is.null(track@size)) "default" else track@size))
            if (!is.null(track_data)) {
                message(sprintf("    Data range: %.2f - %.2f",
                              min(track_data$score),
                              max(track_data$score)))
            }
        }

        list(
            success = TRUE,
            data = track,
            error = NULL
        )

    }, error = function(e) {
        if (verbose) {
            message("  ERROR creating control track:")
            message(sprintf("    %s", e$message))
        }

        list(
            success = FALSE,
            data = NULL,
            error = as.character(e)
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
    #ylim,
    title = NULL,
    visualization_params = list(),
    verbose = FALSE
) {
    # Input validation
    stopifnot(
        "tracks must be a list" = is.list(tracks),
        "tracks cannot be empty" = length(tracks) > 0,
        "chromosome must be character" = is.character(chromosome),
        "from must be numeric" = is.numeric(from),
        "to must be numeric" = is.numeric(to),
       ## "ylim must be numeric vector of length 2" =
       ##     is.numeric(ylim) && length(ylim) == 2,
        "visualization_params must be a list" = is.list(visualization_params)
    )

    if (verbose) {
        message("\nCreating Plot Configuration:")
        message(sprintf("  Tracks: %d", length(tracks)))
        message(sprintf("  Chromosome: %s", chromosome))
        message(sprintf("  Range: %d-%d", from, to))
        #message(sprintf("  Y-limits: %s", paste(ylim, collapse="-")))
    }

    # Validate tracks
    valid_tracks <- all(sapply(tracks, inherits, "GdObject"))
    if (!valid_tracks) {
        stop("All elements in 'tracks' must be valid Gviz track objects")
    }

    # Create base configuration
    base_config <- list(
        trackList = tracks,
        chromosome = chromosome,
        from = from,
        to = to,
        #ylim = ylim,
        main = title
    )

    # Merge with visualization parameters
    final_config <- modifyList(base_config, visualization_params)

    if (verbose) {
        message("\nConfiguration Parameters:")
        message("  Base parameters:")
        invisible(lapply(names(base_config), function(param) {
            message(sprintf("    %s", param))
        }))
        message("  Visualization parameters:")
        invisible(lapply(names(visualization_params), function(param) {
            message(sprintf("    %s", param))
        }))
    }

    final_config
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
    plot_params = list(),
    display_plot = FALSE,
    verbose = TRUE,
    output_format = "svg"
) {
    # Validate inputs
    stopifnot(
        "plot_config must be a list" = is.list(plot_config),
        "plot_config must contain trackList" = !is.null(plot_config$trackList),
        "save_params must be a list" = is.list(save_params),
        "plot_params must be a list" = is.list(plot_params),
        "output_format must be one of 'svg', 'pdf', or 'png'" = output_format %in% c("svg", "pdf", "png")
    )

    # Check for svglite availability
    svglite_available <- requireNamespace("svglite", quietly = TRUE)

    # Check for rsvg and magick availability
    rsvg_available <- requireNamespace("rsvg", quietly = TRUE)
    magick_available <- requireNamespace("magick", quietly = TRUE)

    # Debug output
    if (verbose) {
        message("\nPlot Configuration:")
        message(sprintf("  Number of tracks: %d", length(plot_config$trackList)))
        message("  Track types:")
        invisible(lapply(plot_config$trackList, function(track) {
            message(sprintf("    - %s: %s",
                          class(track)[1],
                          track@name))
        }))

        if (length(plot_params) > 0) {
            message("\nPlot Parameters:")
            invisible(lapply(names(plot_params), function(param) {
                message(sprintf("    %s: %s",
                              param,
                              paste(plot_params[[param]], collapse = ", ")))
            }))
        }

        message("\nDevice Information:")
        message(sprintf("  svglite available: %s", svglite_available))
        message(sprintf("  rsvg available: %s", rsvg_available))
        message(sprintf("  magick available: %s", magick_available))
        message(sprintf("  Output format: %s", output_format))
        message(sprintf("  Display Plot: %s", display_plot))
        message(sprintf("  Save Plot: %s", !is.null(save_path)))
    }

    # Create plotting function
    do_plot <- function() {
        if (verbose) message("\nExecuting plot...")

        final_params <- modifyList(plot_config, plot_params)

        tryCatch({
            do.call(Gviz::plotTracks, final_params)
            if (verbose) message("  Plot execution successful")
        }, error = function(e) {
            message("  Plot execution failed:")
            message(sprintf("    %s", as.character(e)))
            stop(e)
        })
    }

    # Handle saving
    if (!is.null(save_path)) {
        tryCatch({
            # Validate save directory
            save_dir <- dirname(save_path)
            if (!dir.exists(save_dir)) {
                stop("Directory does not exist: ", save_dir)
            }
            if (file.access(save_dir, mode = 2) != 0) {
                stop("Directory is not writable: ", save_dir)
            }
            save_path <- normalizePath(save_path, mustWork = FALSE)

            # Set default save parameters
            default_save_params <- list(
                width = 10,
                height = 8,
                device = ifelse(svglite_available, "svglite", "svg")
            )
            save_params <- modifyList(default_save_params, save_params)

            if (verbose) {
                message("\nSaving plot:")
                message(sprintf("  Path: %s", save_path))
                message(sprintf("  Width: %d", save_params$width))
                message(sprintf("  Height: %d", save_params$height))
                message(sprintf("  Device: %s", save_params$device))
            }

            # Handle different file formats
            if (output_format == "svg") {
                if (svglite_available && save_params$device == "svglite") {
                    if (verbose) message("  Using svglite package for SVG output")
                    svglite::svglite(save_path,
                                    width = save_params$width,
                                    height = save_params$height)
                } else {
                    if (verbose) message("  Using base SVG device")
                    svg(save_path,
                       width = save_params$width,
                       height = save_params$height)
                }
                do_plot()
                dev.off()
            } else if (output_format == "pdf") {
                pdf(save_path,
                    width = save_params$width,
                    height = save_params$height)
                do_plot()
                dev.off()
            } else if (output_format == "png") {
                # Use rsvg or magick to convert SVG to PNG
                if (rsvg_available || magick_available) {
                    # Save as SVG first
                    temp_svg <- tempfile(fileext = ".svg")
                    if (svglite_available && save_params$device == "svglite") {
                        svglite::svglite(temp_svg,
                                        width = save_params$width,
                                        height = save_params$height)
                    } else {
                        svg(temp_svg,
                           width = save_params$width,
                           height = save_params$height)
                    }
                    do_plot()
                    dev.off()

                    # Convert SVG to PNG
                    if (rsvg_available) {
                        if (verbose) message("  Using rsvg to convert SVG to PNG")
                        rsvg::rsvg_png(temp_svg, save_path, width = save_params$width * 100, height = save_params$height * 100)
                    } else if (magick_available) {
                        if (verbose) message("  Using magick to convert SVG to PNG")
                        magick::image_write(magick::image_read_svg(temp_svg), save_path, format = "png")
                    }
                    file.remove(temp_svg)  # Clean up temporary SVG file
                } else {
                    stop("Neither rsvg nor magick is available for PNG conversion. Please install one of them.")
                }
            } else {
                stop("Unsupported output format. Use 'svg', 'pdf', or 'png'")
            }
        }, error = function(e) {
            stop("Failed to save plot: ", e$message)
        })
    }

    # Handle display
    if (display_plot) {
        if (verbose) message("\nDisplaying plot...")

        if (dev.cur() == 1) {
            if (verbose) message("  Opening new graphics device")
            if (capabilities("aqua")) {
                quartz()
            } else if (.Platform$OS.type == "windows") {
                windows()
            } else {
                X11()
            }
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
