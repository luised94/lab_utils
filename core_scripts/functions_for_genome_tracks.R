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
#' @description Creates an empty DataTrack with evenly spaced zero values across
#'   a specified chromosome range
#'
#' @param sampling_rate numeric Number of base pairs per data point
#' @param chromosome_width numeric Total width of chromosome in base pairs
#' @param track_color character Color for track visualization
#' @param type character Track type ('l' for line, 'h' for histogram, etc.)
#' @param chromosome_name character Chromosome identifier
#' @param placeholder_format_name character Format string for track name
#' @param ... Additional arguments passed to sprintf for track naming
#'
#' @return Gviz DataTrack object with placeholder data
#'
#' @examples
#' create_placeholder_track(
#'   sampling_rate = 100,
#'   chromosome_width = 1000,
#'   track_color = "#CCCCCC",
#'   chromosome_name = "chrI",
#'   placeholder_format_name = "%s (No data)"
#' )
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom Gviz DataTrack
create_placeholder_track <- function(
    sampling_rate = 100,
    chromosome_width = NULL,
    track_color,
    track_type = "l",
    chromosome_name,
    placeholder_format_name,
    ...
) {
    # Input validation
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
            is.character(placeholder_format_name)
    )

    # Create placeholder track (your existing logic)
    num_points <- ceiling(chromosome_width / sampling_rate)
    placeholder_name <- sprintf(placeholder_format_name, ...)

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
        type = track_type,
        col = track_color,
        chromosome = chromosome_name
    )
}

create_sample_track <- function(
    bigwig_file_path,
    track_format_name,
    sample_track_name,
    track_color,
    track_type,
    genomic_range,
    ...
) {
    stopifnot()
    if(!is.na(bigwig_file_path) && file.exists(bigwig_file)){

    track_data <- rtracklayer::import(bigwig_file, which = genomic_range)
    track_name <- sprintf(track_format_name, ...)

    return(Gviz::DataTrack(
        track_data,
        name = track_name,
        type = "l",
        col = track_color
    ))
    } else {
        # Not sure if it should be logical, null, string, etcetera.
        return(NULL)
    }

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
