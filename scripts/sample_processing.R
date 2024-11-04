# Validate and prepare categories for labeling
sample_validate_categories <- function(table, categories) {
    if (!is.data.frame(table)) {
        stop("Input table must be a data frame")
    }
    if (!is.character(categories) || length(categories) == 0) {
        stop("Categories must be a non-empty character vector")
    }
    
    # Ensure antibody category is included
    if (!"antibody" %in% categories) {
        categories <- c("antibody", categories)
    }
    
    # Validate categories exist in table
    missing_categories <- setdiff(categories, colnames(table))
    if (length(missing_categories) > 0) {
        stop("Missing categories in table: ", 
             paste(missing_categories, collapse = ", "))
    }
    
    return(categories)
}

# Get unique values for each category
sample_get_unique_values <- function(table, categories) {
    unique_values <- lapply(table[categories], unique)
    names(unique_values) <- categories
    return(unique_values)
}

# Create label for a single sample
sample_create_label <- function(sample, categories, unique_values) {
    relevant_categories <- sapply(categories, function(cat) {
        if (length(unique_values[[cat]]) > 1 || cat == "antibody") {
            return(as.character(sample[cat]))
        }
        return(NULL)
    })
    
    # Filter out NULL values and combine
    valid_categories <- relevant_categories[!sapply(relevant_categories, is.null)]
    return(paste(valid_categories, collapse = "_"))
}

# Main function to generate unique labels
sample_generate_labels <- function(table, categories_for_label, verbose = FALSE) {
    # Validate and prepare categories
    validated_categories <- sample_validate_categories(table, categories_for_label)
    
    # Get unique values
    unique_vals <- sample_get_unique_values(table, validated_categories)
    
    # Generate labels for each row
    labels <- apply(table, 1, function(row) {
        sample_create_label(row, validated_categories, unique_vals)
    })
    
    if (verbose) {
        message("Categories used: ", paste(validated_categories, collapse = ", "))
        message("Unique values by category:")
        print(unique_vals)
        message("Generated labels:")
        print(labels)
    }
    
    return(unname(labels))
}
