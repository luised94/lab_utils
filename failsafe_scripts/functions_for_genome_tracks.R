create_minimal_identifiers_alt <- function(sample_ids) {
    # Validation
    stopifnot(
        "sample_ids must be character vector" = is.character(sample_ids),
        "sample_ids cannot be empty" = length(sample_ids) > 0,
        "sample_ids must be unique" = !any(duplicated(sample_ids))
    )
    # Find positions where values differ
    id_matrix <- do.call(rbind, strsplit(sample_ids, ""))
    diff_positions <- which(apply(id_matrix, 2, function(x) length(unique(x)) > 1))
    
    # Get minimal required positions
    min_pos <- min(diff_positions)
    max_pos <- max(diff_positions)
    
    # Extract minimal substring that ensures uniqueness
    short_ids <- substr(sample_ids, min_pos, max_pos + 1)
    
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
