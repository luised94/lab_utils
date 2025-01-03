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
find_control_sample <- function(experimental_sample, metadata, control_factors, verbose) {
    # Create matching conditions for control
    control_conditions <- lapply(control_factors$genotype, function(factor) {
        metadata[[factor]] == experimental_sample[[factor]]
    })
    
    # Combine conditions with Input antibody requirement
    control_conditions$is_input <- metadata$antibody == "Input"
    
    # Find matching control samples
    control_matches <- Reduce(`&`, control_conditions)
    
    if (sum(control_matches) == 0) {
        if (verbose) {
            message("No matching control found for sample: ", 
                   experimental_sample$sample_id)
        }
        return(NULL)
    }
    
    # Return first matching control
    metadata[control_matches, ][1, ]
}

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

#' @title Create Track Labels from Sample Metadata
#' @description Creates labels based on sample categories with configurable display rules
#' @param samples data.frame Sample metadata
#' @param categories character vector Categories to consider for labels
#' @param always_show character vector Categories to always show
#' @param never_show character vector Categories to never show
#' @param separator character String to use between label components
#' @param verbose logical Print processing information
create_track_labels <- function(samples, 
                              categories = NULL,  # Now optional
                              always_show = c("antibody"),
                              never_show = c("sample_id", "full_name", "short_name", "X__cf_genotype"),
                              separator = " - ",
                              verbose = FALSE) {
    result <- tryCatch({
        # Input validation
        stopifnot(
            "samples must be data.frame" = is.data.frame(samples),
            "always_show must be character" = is.character(always_show),
            "never_show must be character" = is.character(never_show),
            "separator must be character" = is.character(separator),
            "always_show categories must exist in samples" = 
                all(always_show %in% colnames(samples))
        )
        
        # If categories not provided, use all columns except never_show
        if (is.null(categories)) {
            categories <- setdiff(colnames(samples), never_show)
            if (verbose) {
                message("Using all available categories except: ", 
                       paste(never_show, collapse = ", "))
            }
        } else {
            stopifnot("categories must be character" = is.character(categories))
            # Remove any never_show categories if present
            categories <- setdiff(categories, never_show)
        }
        
        if (verbose) {
            message("Analyzing categories: ", paste(categories, collapse = ", "))
        }
        
        # Find categories with varying values
        varying_categories <- categories[sapply(categories, function(cat) {
            length(unique(samples[[cat]])) > 1
        })]
        
        if (verbose) {
            message("Found varying categories: ", 
                   paste(varying_categories, collapse = ", "))
        }
        
        # Combine always_show with varying categories, ensuring order
        label_categories <- unique(c(
            always_show,  # Always first
            intersect(varying_categories, categories)  # Then varying
        ))
        
        # Remove any never_show categories that might have slipped through
        label_categories <- setdiff(label_categories, never_show)
        
        if (verbose) {
            message("Final label categories: ", 
                   paste(label_categories, collapse = ", "))
        }
        
        # Create category information
        category_info <- list(
            distinguishing = label_categories,
            varying = varying_categories,
            always_shown = always_show,
            never_shown = never_show,
            all_available = colnames(samples),
            used_separator = separator
        )
        
        # Create labels
        labels <- apply(samples, 1, function(row) {
            # Get values for each category
            values <- sapply(label_categories, function(cat) {
                value <- row[[cat]]
                # Clean up NA values
                if (is.na(value)) return("NA")
                as.character(value)
            })
            
            # Remove empty or NA values
            values <- values[values != "" & values != "NA"]
            
            # Combine with separator
            paste(values, collapse = separator)
        })
        
        if (verbose) {
            message("\nLabel Summary:")
            message("- Total labels: ", length(labels))
            message("- Unique labels: ", length(unique(labels)))
            message("- Example label: ", labels[1])
        }
        
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
