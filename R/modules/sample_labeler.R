#' Sample Labeling Functions
unique_labeling <- function(table,
                          categories_for_label,
                          config = CONFIG$LABEL_CONFIG) {
    log_info("Creating unique labels for samples")
    
    # Validate inputs
    validate_labeling_inputs(table, categories_for_label)
    
    # Ensure required categories
    categories_for_label <- ensure_required_categories(categories_for_label)
    
    # Get unique values
    unique_values <- get_unique_category_values(table, categories_for_label)
    
    # Create labels
    labels <- create_sample_labels(table,
                                 categories_for_label,
                                 unique_values,
                                 config)
    
    return(labels)
}

validate_labeling_inputs <- function(table, categories) {
    if (!is.data.frame(table)) {
        log_error("Invalid input table type")
        stop("Table must be a data frame")
    }
    
    if (!is.character(categories) || length(categories) == 0) {
        log_error("Invalid categories specification")
        stop("Categories must be non-empty character vector")
    }
    
    missing_cats <- setdiff(categories, colnames(table))
    if (length(missing_cats) > 0) {
        log_error("Missing categories:", paste(missing_cats, collapse = ", "))
        stop("Missing categories in table")
    }
}

ensure_required_categories <- function(categories) {
    required <- CONFIG$REQUIRED_CATEGORIES$BASE
    if (!required %in% categories) {
        log_info("Adding required category:", required)
        categories <- c(required, categories)
    }
    return(categories)
}

create_sample_labels <- function(table,
                               categories,
                               unique_values,
                               config) {
    log_info("Constructing sample labels")
    
    labels <- apply(table, 1, function(sample) {
        # Get relevant category values
        values <- sapply(categories, function(cat) {
            if (length(unique_values[[cat]]) > 1 || 
                cat == CONFIG$REQUIRED_CATEGORIES$BASE) {
                return(sample[cat])
            }
            return(NULL)
        })
        
        # Filter and combine values
        values <- values[!sapply(values, is.null)]
        label <- paste(values, collapse = config$SEPARATOR)
        
        # Truncate if necessary
        if (nchar(label) > config$MAX_LENGTH) {
            label <- paste0(substr(label, 1,
                                 config$MAX_LENGTH - nchar(config$TRUNCATE_SUFFIX)),
                          config$TRUNCATE_SUFFIX)
        }
        
        return(label)
    })
    
    return(unlist(labels))
}
