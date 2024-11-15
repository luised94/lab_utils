create_minimal_identifiers <- function(sample_ids, verbose = FALSE) {

    # Validation
    stopifnot(
        "sample_ids must be character vector" = is.character(sample_ids),
        "sample_ids cannot be empty" = length(sample_ids) > 0,
        "sample_ids must be unique" = !any(duplicated(sample_ids)),
        "sample_ids must have equal length" = length(unique(nchar(sample_ids))) == 1,
        "verbose must be logical" = is.logical(verbose)
    )

    if (verbose) {
        message(sprintf("Processing %d sample IDs of length %d", 
                       length(sample_ids), nchar(sample_ids[1])))
    }
    # Find positions where values differ
    id_matrix <- do.call(rbind, strsplit(sample_ids, ""))
    diff_positions <- which(apply(id_matrix, 2, function(x) length(unique(x)) > 1))

    if (verbose) {
        message(sprintf("Found differences at positions: %s", 
                       paste(diff_positions, collapse = ", ")))
    }
    
    # Get minimal required positions
    min_pos <- min(diff_positions)
    max_pos <- max(diff_positions)
    
    # Extract minimal substring that ensures uniqueness
    short_ids <- substr(sample_ids, min_pos, max_pos + 1)

    if (verbose) {
        message(sprintf("Reduced sample IDs from %d to %d digits", 
                       nchar(sample_ids[1]), nchar(short_ids[1])))
        message(sprintf("Using positions %d to %d", min_pos, max_pos + 1))
    }
    
    # Verify uniqueness of result
    if (any(duplicated(short_ids))) {
        stop("Failed to create unique short identifiers")
    }
    
    return(short_ids)
}

generate_distinct_colors <- function(n) {
    # RColorBrewer provides good distinct colors
    if (n <= 8) {
        RColorBrewer::brewer.pal(max(3, n), "Set2")[1:n]
    } else {
        # For more categories, use rainbow with better spacing
        rainbow(n, s = 0.7, v = 0.9)
    }
}

# Function to find matching control sample
find_control_sample <- function(experimental_sample, metadata, control_factors) {
    # Create matching conditions for control
    control_conditions <- lapply(control_factors$genotype, function(factor) {
        metadata[[factor]] == experimental_sample[[factor]]
    })
    
    # Combine conditions with Input antibody requirement
    control_conditions$is_input <- metadata$antibody == "Input"
    
    # Find matching control samples
    control_matches <- Reduce(`&`, control_conditions)
    
    if (sum(control_matches) == 0) {
        if (DEBUG_CONFIG$verbose) {
            message("No matching control found for sample: ", 
                   experimental_sample$sample_id)
        }
        return(NULL)
    }
    
    # Return first matching control
    metadata[control_matches, ][1, ]
}

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


#' @title Create Track Labels and Category Information
#' @description Creates labels and returns category usage information
#' @param samples data.frame Sample metadata
#' @param categories character vector Categories to consider
#' @param always_show character vector Categories to always show
#' @param verbose logical Print processing information
#' @return list containing labels and category information
create_track_labels <- function(samples, categories, 
                              always_show = c("antibody"), 
                              verbose = FALSE) {
    result <- tryCatch({
        # Input validation
        stopifnot(
            "samples must be data.frame" = is.data.frame(samples),
            "categories must be character" = is.character(categories),
            "categories cannot be empty" = length(categories) > 0,
            "always_show must be character" = is.character(always_show),
            "categories must exist in samples" = 
                all(c(categories, always_show) %in% colnames(samples))
        )
        
        if (verbose) {
            message("Creating labels using categories: ", 
                   paste(categories, collapse = ", "))
        }
        
        # Find distinguishing categories
        varying_categories <- categories[sapply(categories, function(cat) {
            length(unique(samples[[cat]])) > 1
        })]
        
        # Combine with always_show categories
        distinguishing_cats <- unique(c(always_show, varying_categories))
        
        # Track category usage
        category_info <- list(
            distinguishing = distinguishing_cats,
            varying = varying_categories,
            always_shown = always_show,
            all_available = colnames(samples)
        )
        
        if (verbose) {
            message("Distinguishing categories: ", 
                   paste(distinguishing_cats, collapse = ", "))
        }
        
        # Create labels
        labels <- apply(samples, 1, function(row) {
            values <- sapply(distinguishing_cats, function(cat) row[[cat]])
            paste(values, collapse = " - ")
        })
        
        list(
            success = TRUE,
            data = list(
                labels = labels,
                categories = category_info
            ),
            error = NULL
        )
    }, error = function(e) {
        list(
            success = FALSE,
            data = NULL,
            error = e$message
        )
    })
    
    return(result)
}

#' @title Create Consistent Color Scheme
#' @description Generates and manages consistent colors for track visualization
#' @param config list Configuration including fixed colors
#' @param categories list Named list of category values requiring colors
#' @return list Color assignments and management functions
create_color_scheme <- function(config, categories, verbose = FALSE) {
    # Input validation
    stopifnot(
        "config must be a list" = is.list(config),
        "categories must be a list" = is.list(categories),
        "config must contain placeholder color" = !is.null(config$placeholder),
        "config must contain input color" = !is.null(config$input),
        "verbose must be logical" = is.logical(verbose)
    )
    
    # Initialize fixed colors
    fixed_colors <- list(
        placeholder = config$placeholder,
        input = config$input
    )
    
    if (verbose) {
        message("Initializing color scheme...")
        message("Fixed colors:")
        message("- Placeholder: ", fixed_colors$placeholder)
        message("- Input: ", fixed_colors$input)
    }
    
    # Generate category-specific colors
    category_colors <- list()
    
    for (category_name in names(categories)) {
        # Get unique values for current category
        category_values <- categories[[category_name]]
        color_count <- length(category_values)
        
        # Set consistent seed for reproducibility
        category_seed <- sum(utf8ToInt(category_name))
        set.seed(category_seed)
        
        # Generate colors based on category size
        if (color_count <= 8) {
            category_palette <- RColorBrewer::brewer.pal(max(3, color_count), "Set2")
            category_colors_vector <- category_palette[seq_len(color_count)]
        } else {
            category_colors_vector <- rainbow(color_count, s = 0.7, v = 0.9)
        }
        
        # Assign names to colors
        names(category_colors_vector) <- category_values
        category_colors[[category_name]] <- category_colors_vector
        
        if (verbose) {
            message(sprintf("\nColors for category '%s':", category_name))
            for (value in category_values) {
                message(sprintf("- %s: %s", value, category_colors[[category_name]][value]))
            }
        }
    }
    
    # Create color getter function
    get_color <- function(category, value) {
        if (category == "antibody" && value == "Input") {
            return(fixed_colors$input)
        }
        
        if (!category %in% names(category_colors)) {
            if (verbose) {
                message(sprintf("Category '%s' not found, using placeholder", category))
            }
            return(fixed_colors$placeholder)
        }
        
        if (!value %in% names(category_colors[[category]])) {
            if (verbose) {
                message(sprintf("Value '%s' not found in category '%s', using placeholder",
                              value, category))
            }
            return(fixed_colors$placeholder)
        }
        
        return(category_colors[[category]][value])
    }
    
    # Create output structure
    color_scheme <- list(
        fixed = fixed_colors,
        categories = category_colors,
        get_color = get_color
    )
    
    return(color_scheme)
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
        col1 <- sprintf(
            "Experiment: %s\nComparison: %s\nChromosome: %s\nSamples: %d\nNorm: %s",
            plot_info$experiment_id,
            comparison_name,
            plot_info$chromosome,
            nrow(metadata),
            plot_info$normalization
        )
        
        col2 <- if (length(shared_values) > 0) {
            sprintf("Shared Properties:\n%s", paste(shared_values, collapse = "\n"))
        } else {
            "No Shared Properties"
        }
        
        col3 <- sprintf(
            "Processing:\n%s\nY-range: [%.2f, %.2f]",
            format(Sys.time(), "%Y-%m-%d %H:%M"),
            plot_info$y_limits[1],
            plot_info$y_limits[2]
        )
        
        title <- sprintf("%-40s %-40s %-40s", col1, col2, col3)
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
