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
    format_args
) {
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

    num_points <- ceiling(chromosome_width / sampling_rate)
    placeholder_name <- do.call(sprintf, c(list(placeholder_format_name), format_args))

    empty_ranges <- GenomicRanges::GRanges(
        seqnames = chromosome_name,
        ranges = IRanges::IRanges(
            start = seq(1, chromosome_width, length.out = num_points),
            width = 1
        ),
        score = rep(0, num_points),
        strand = "*"
    )
    
    Gviz::DataTrack(
        empty_ranges,
        name = placeholder_name,
        type = type,
        col = track_color,
        chromosome = chromosome_name
    )
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
    genomic_range = NULL
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

    # Check file existence
    if (is.na(bigwig_file_path) || !file.exists(bigwig_file_path)) {
        return(list(
            success = FALSE,
            data = NULL,
            error = "Bigwig file not found"
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
