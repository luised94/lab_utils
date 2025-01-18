#' Print Configuration Settings
#' @param config List containing configuration settings
#' @param title Optional custom title for the configuration display
#' @param verbose Boolean to control whether to print
#' @return Invisible return of the formatted settings
print_config_settings <- function(config, title = NULL, verbose = TRUE) {
    # Input validation
    stopifnot(
        "config must be a list" = is.list(config),
        "config cannot be empty" = length(config) > 0,
        "verbose must be logical" = is.logical(verbose) && length(verbose) == 1
    )

    if (!verbose) {
        return(invisible(config))
    }

    # Format title
    if (is.null(title)) {
        title <- deparse(substitute(config))
    }
    
    # Create formatted output
    header <- sprintf("=== %s Settings ===", title)
    separator <- paste(rep("=", nchar(header)), collapse = "")
    
    cat(sprintf("\n%s\n", header))
    
    # Process and format each setting
    settings <- lapply(names(config), function(setting) {
        value <- config[[setting]]
        
        # Format different types of values
        formatted_value <- if (is.logical(value)) {
            if(value) "YES" else "NO"
        } else if (is.null(value)) {
            "NULL"
        } else if (is.list(value)) {
            "LIST"  # Could expand this for nested configs if needed
        } else if (is.vector(value) & length(value) > 0) {
            # Handle vectors (like column order) specially
            paste0("[", paste(value, collapse = ", "), "]")
        } else {
            stop(sprintf(
                "Unsupported configuration value type for setting '%s':\n  Type: %s\n  Value: %s",
                setting,
                class(value)[1],
                deparse(value)
            ))
        }
        
        # Format the line
        sprintf("%-25s: %s", 
                gsub("_", " ", toupper(setting)), 
                formatted_value)
    })
    
    # Print settings
    cat(paste(settings, collapse = "\n"), "\n")
    cat(separator, "\n")
    
    # Return invisibly
    invisible(config)
}
