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
        max_width = PLOT_CONFIG$title$format$max_width,
        max_lines = PLOT_CONFIG$title$format$max_lines
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
## Update PLOT_CONFIG with dynamic colors
#PLOT_CONFIG$track_colors <- list(
#    antibody = antibody_colors,
#    placeholder = PLOT_CONFIG$placeholder_color  # Maintain consistent placeholder
#)

# Create color legend text
#legend_text <- sprintf(
#    "Track Colors:\n%s\n%s",
#    paste("?", names(PLOT_CONFIG$track_colors$antibody), 
#          sprintf("(%s)", PLOT_CONFIG$track_colors$antibody), 
#          collapse = "\n"),
#    sprintf("? No Data (%s)", PLOT_CONFIG$placeholder_color)
#)
