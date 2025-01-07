validate_configs <- function(required_configs) {
    missing <- required_configs[!sapply(required_configs, exists)]
    if (length(missing) > 0) {
        stop("Missing required configs: ", paste(missing, collapse=", "))
    }
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
