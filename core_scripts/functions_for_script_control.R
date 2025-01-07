library(optparse)

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

#' Parse common command line arguments for experiment scripts
#' @param description Character string describing the script's purpose
#' @param prog Character string of program name (optional)
#' @return List of validated arguments
#' @examples
#' args <- parse_common_arguments("Align sequencing reads to reference genome")
parse_common_arguments <- function(description = "Script description", prog = NULL) {
    option_list <- list(
        make_option(
            "--experiment-id",
            type = "character",
            help = "Experiment ID (format: YYMMDD'Bel', e.g., 241228Bel)",
            dest = "experiment_id",
            metavar = "ID"
        ),
        make_option(
            "--confirmation-only",
            action = "store_true",
            default = TRUE,
            help = "Only validate and display configuration, then exit",
            dest = "confirmation_only"
        )
    )

    # Create parser with optional program name
    opt_parser <- OptionParser(
        option_list = option_list,
        description = description,
        prog = prog
    )

    # Parse arguments with enhanced error handling
    args <- tryCatch(
        parse_args(opt_parser),
        error = function(e) {
            cat("\nERROR: Argument parsing failed\n", file = stderr())
            cat(as.character(e), "\n\n", file = stderr())
            print_help(opt_parser)
            quit(status = 1, save = "no")
        }
    )

    # Validate required argument
    if (is.null(args$experiment_id)) {
        print_help(opt_parser)
        stop("Missing required argument: --experiment-id")
    }

    # Validate experiment ID format
    if (!grepl("^\\d{6}Bel$", args$experiment_id, perl = TRUE)) {
        stop(sprintf(
            "Invalid experiment-id format.\nExpected: YYMMDD'Bel'\nReceived: %s",
            args$experiment_id
        ))
    }

    return(args)
}
