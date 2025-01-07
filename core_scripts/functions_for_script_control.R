validate_configs <- function(required_configs) {
    # 1. Input Validation
    stopifnot(
        "required_configs must be a character vector" = 
            is.character(required_configs),
        "required_configs cannot be empty" = 
            length(required_configs) > 0,
        "required_configs cannot contain NA values" = 
            !any(is.na(required_configs)),
        "required_configs cannot contain empty strings" = 
            !any(required_configs == ""),
        "required_configs must be unique" = 
            !any(duplicated(required_configs))
    )

    # 2. Ensure each item ends with "_CONFIG"
    invalid_pattern <- grep(".*_CONFIG$", required_configs, invert = TRUE, value = TRUE)
    if (length(invalid_pattern) > 0) {
        stop(
          "All required_configs must end with '_CONFIG'. Invalid entries: ",
          paste(invalid_pattern, collapse = ", ")
        )
    }
    
    # 3. Verify each one actually exists in the environment
    missing <- required_configs[!sapply(required_configs, exists)]
    if (length(missing) > 0) {
        stop("Missing required configs: ", paste(missing, collapse = ", "))
    }
    
    # If we get here, all checks have passed
    invisible(TRUE)
}

#
parse_args <- function(args) {
    # Validate args input characteristics
    stopifnot(
        "args must be a character vector" = is.character(args),
        "args must come from commandArgs" = {
            # Check typical commandArgs characteristics
            is.atomic(args) && 
            is.vector(args) && 
            !is.null(attr(args, "names")) == FALSE && # Should not have names
            length(attributes(args)) <= 1  # Should have minimal attributes
        }
    )
    
    # Handle empty args
    if (length(args) == 0) {
        stop(
            "No arguments provided\n",
            "Usage: Rscript script.R --experiment-id=YYMMDD'Bel'\n",
            "Example: --experiment-id=241228Bel"
        )
    }
    
    # Rest of the function remains the same
    args_chunks <- strsplit(paste(args, collapse=" "), " ")[[1]]
    args_list <- list()
    
    for (arg in args_chunks) {
        # Handle non-value flags
        if (!grepl("=", arg)) {
            key <- sub("^--", "", arg)
            args_list[[key]] <- TRUE
            next
        }
        
        # Handle key=value pairs
        parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
        if (length(parts) != 2) {
            stop("Invalid argument format: ", arg, "\nExpected format: --key=value")
        }
        key <- sub("^--", "", parts[1])
        value <- parts[2]
        
        # Smart type conversion
        value <- if (tolower(value) %in% c("true", "false")) {
            tolower(value) == "true"
        } else if (grepl("^[0-9]+$", value)) {
            as.integer(value)
        } else if (grepl("^[0-9]*\\.?[0-9]+$", value)) {
            as.numeric(value)
        } else {
            value
        }
        
        args_list[[key]] <- value
    }
    
    # Validate experiment-id
    if (is.null(args_list[["experiment-id"]])) {
        stop("Missing required argument: --experiment-id")
    }
    
    if (!grepl("^\\d{6}Bel$", args_list[["experiment-id"]], perl = TRUE)) {
        stop("Invalid experiment-id format. Must be 6 digits followed by 'Bel', got: ", 
             args_list[["experiment-id"]])
    }
    
    args_list
}

#' Prompt for interactive configuration confirmation
#' @return Logical indicating whether to proceed
#' @examples
#' if (!confirm_configuration()) {
#'     stop("Script terminated by user")
#' }
confirm_configuration <- function() {
    if (!interactive()) {
        stop("This script requires interactive confirmation")
    }
    
    cat("\nPlease review the configuration settings above carefully.\n")
    
    while (TRUE) {
        response <- readline("Proceed with these settings? [y/N]: ")
        
        # Default to No on empty input
        if (response == "") {
            message("Configuration rejected (default: No)")
            return(FALSE)
        }
        
        # Check response
        if (tolower(response) %in% c("y", "yes")) {
            return(TRUE)
        } else if (tolower(response) %in% c("n", "no")) {
            message("Configuration rejected by user")
            return(FALSE)
        }
        
        message("Invalid input. Please answer 'yes' or 'no' (or 'y' or 'n')")
    }
}
