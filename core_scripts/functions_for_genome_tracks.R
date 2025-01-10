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

